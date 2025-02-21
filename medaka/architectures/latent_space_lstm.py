"""Calculates read level embeddings, then pools by averaging embeddings."""
import warnings

import torch

import medaka.architectures.base_classes as base_classes
import medaka.architectures.read_level_modules as read_level_modules
import medaka.features


class ReversibleLSTM(torch.nn.Module):
    """A wrapper around torch.nn.LSTM that allows for reversed layers."""

    def __init__(self, *args, reverse=False, **kwargs):
        """Initialise the reversible LSTM."""
        super().__init__()
        self.lstm = torch.nn.LSTM(*args, **kwargs)
        self.reverse = reverse

    def forward(self, x):
        """Forward pass of the LSTM."""
        flip_dim = 1 if self.lstm.batch_first else 0
        if self.reverse:
            x = x.flip(flip_dim)
        x = self.lstm.forward(x)[0]
        if self.reverse:
            x = x.flip(flip_dim)
        return x

    def __repr__(self):
        """Return string representation."""
        s = super().__repr__()
        return s + f", reversed={self.reverse})"


class LatentSpaceLSTM(base_classes.ReadLevelFeaturesModel):
    """Implementation of latent space LSTM model.

    Takes read level features (see medaka.features.ReadAlignmentFeatureEncoder)
    as input and maps onto haploid per position base predictions.

    NOTE: Currently this model uses only basecalls, quality scores, strand and
    (if use_dwells=True) the dwells channel. This means that the mapQC, channel
    and dtype options are all ignored.
    """

    def __init__(
        self,
        num_classes=5,
        lstm_size=128,
        cnn_size=128,
        kernel_sizes=[1, 17],
        pooler_type="mean",
        pooler_args={},
        use_dwells=False,
        bases_alphabet_size=6,
        bases_embedding_size=6,
        bidirectional=True,
        time_steps=None,
    ):
        """Initialise bidirectional latent space LSTM model.

        :param num_classes: int, number of output class labels.
        :param lstm_size: int, size of each lstm layer (in each direction).
        :param cnn_size: int, number of channels in the read conv layers.
        :param kernel_sizes: List[int], kernel sizes for read level conv layers
        :param pooler_type: str, type of layer to pool read level embeddings
            (see read_level_modules.Poolers for options)
        :param pooler_args: dict, arguments for the pooling layer
        :param use_dwells: bool, if True, use dwell times as an extra feature
        :param bases_alphabet_size: int, number of bases in the alphabet
        :param bases_embedding_size: int, size of the base embeddings
        :param bidirectional: bool, if True, use 2-layer bidirectional LSTM,
            else 4 interleaved unidirectional layers are used
        :param time_steps: dummy argument, no longer used
        """
        super().__init__()

        if time_steps is not None:
            warnings.warn("timesteps is no lnoger required to be specified")

        # set properties needed for export to dictionary,
        # see TorchModel.to_dict
        self.num_classes = num_classes
        self.lstm_size = lstm_size
        self.cnn_size = cnn_size
        self.kernel_sizes = kernel_sizes
        self.pooler_type = pooler_type
        self.pooler_args = pooler_args
        self.use_dwells = use_dwells
        self.bases_alphabet_size = bases_alphabet_size
        self.bases_embedding_size = bases_embedding_size

        # embedding for bases
        self.base_embedder = torch.nn.Embedding(
            bases_alphabet_size, bases_embedding_size
        )
        # strand can be [backward (-1), forward (1), or none (0)]
        self.strand_embedder = torch.nn.Embedding(3, bases_embedding_size)

        # find the number of extra features that will be concatenated onto
        # the bases and strand embeddings
        extra_features = 1  # quality scores
        if self.use_dwells:
            extra_features += 1

        # make read level convolutional layers
        self.read_level_conv = read_level_modules.ReadLevelConv(
            in_features=bases_embedding_size + extra_features,
            out_dim=self.lstm_size,
            kernel_sizes=self.kernel_sizes,
            channel_dim=self.cnn_size,
            use_batch_norm=True,
        )

        # linear projection layer to expand the output of the read level conv
        # to be the same size as the LSTM input
        self.pre_pool_expansion_layer = torch.nn.Linear(
            self.cnn_size, self.lstm_size
        )

        # make pooling layer
        self.pooler = read_level_modules.Poolers.from_str(self.pooler_type)(
            **self.pooler_args
        )

        # bidirectional LSTM
        self.bidirectional = bidirectional
        if bidirectional:
            self.lstm = torch.nn.LSTM(
                lstm_size,
                lstm_size,
                num_layers=2,
                bidirectional=True,
                batch_first=True,
            )
        # otherwise unidirectional
        else:
            self.lstm = torch.nn.Sequential(
                *[ReversibleLSTM(
                    lstm_size,
                    lstm_size,
                    batch_first=True,
                    reverse=not bool(i % 2)  # reverse-forward-reverse-forward
                ) for i in range(4)]
            )

        # now project the LSTM output to the number of classes
        self.linear = torch.nn.Linear(
            (1 + bidirectional) * lstm_size,
            self.num_classes)
        self.normalise = True

    def forward(self, x):
        """Model forward pass.

        :param x: torch.Tensor, read level feature matrix, shape
                (num_batch, num_positions, num_reads (padded), num_features).

        :return:
            torch.Tensor, logits for positionwise predictions
                (num_positions, num_classes).
        """
        non_empty_position_mask = (
            x.sum((1, -1)) != 0
        )  # shape batch_size x n_reads

        # get the base and strand embeddings for each position in each read
        bases_embedding = self.base_embedder(x[:, :, :, 0].long())
        strand_embedding = self.strand_embedder(x[:, :, :, 2].long() + 1)
        # map the quality scores from approximately [0,50] -> [-1,1]
        # done to ensure roughly equal magnitude of input channels
        scaled_q_scores = (x[:, :, :, 1] / 25 - 1).unsqueeze(-1)
        if self.use_dwells:
            assert (
                x.shape[-1] == 5
            ), "if using dwells, x must have 5 features/read/position"
            dwells = x[:, :, :, 4].unsqueeze(-1)
            x = torch.cat(
                [bases_embedding + strand_embedding, scaled_q_scores, dwells],
                dim=-1,
            )
        else:
            x = torch.cat(
                [bases_embedding + strand_embedding, scaled_q_scores], dim=-1
            )
        # x has shape (batch, positions, reads, features)
        x = x.permute(0, 2, 3, 1)  # batch x reads x features x positions
        b, d, f, p = x.shape
        x = x.flatten(0, 1)
        x = self.read_level_conv(x)  # b*d x cnn_size x p
        x = x.permute(0, 2, 1)  # b*d x p x cnn_size
        x = self.pre_pool_expansion_layer(x)  # b*d x p x lstm_size
        x = x.view(b, d, p, self.lstm_size)
        # mean pool all reads
        x = self.pooler(x, non_empty_position_mask)
        x = self.lstm(x)
        if isinstance(x, tuple):
            x = x[0]  # bi-lstm return (output, hidden)
        x = self.linear(x)

        if self.normalise:
            # switch normalise off for training to return logits for
            # cross-entropy loss
            x = torch.softmax(x, dim=-1)

        return x

    def check_feature_encoder_compatibility(self, fenc):
        """Check feature encoder is valid for this model."""
        clsname = type(self).__name__

        # check is read level feature encoder
        if not isinstance(fenc, medaka.features.ReadAlignmentFeatureEncoder):
            msg = f"{clsname} expects a ReadAlignmentFeatureEncoder."
            raise ValueError(msg)

        # check only one data type
        if len(fenc.dtypes) > 1:
            msg = f"{clsname} is currently only implemented for one dtype."
            raise NotImplementedError(msg)

        if self.use_dwells:
            if not getattr(fenc, "include_dwells", False):
                msg = (
                    "Model expects dwells, however include_dwells not set in "
                    "the feature encoder."
                )
                raise ValueError(msg)

        if getattr(fenc, 'include_haplotypes', False):
            msg = (
                "include_haplotypes attribute set in feature encoder. "
                f"This feature is currently ignored in {clsname}."
            )
            raise self.logger.warn(msg)
