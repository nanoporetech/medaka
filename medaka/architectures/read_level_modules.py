"""Implements models which run 1d convolutions along the read features."""
import enum

import torch


def make_1dconv_layers(
    kernel_sizes, in_feat, channels, use_batch_norm, activation="ReLU"
):
    """Make the 1d conv blocks."""
    if not isinstance(channels, list):
        channels = [channels] * len(kernel_sizes)

    if not isinstance(kernel_sizes, list):
        raise TypeError("kernel_sizes must be a list, one for each conv layer")
    if len(channels) != len(kernel_sizes):
        raise TypeError(
            "if channels is a list, there must be one for each conv layer"
        )

    for k in kernel_sizes:
        msg = "kernel sizes must be odd (for equal & symmetric padding)"
        assert k % 2 == 1, msg

    if activation == "ReLU":
        activation = torch.nn.ReLU
    else:
        raise NotImplementedError(f"Activation {activation} not implemented")

    layers = list()
    for i, (k, c) in enumerate(zip(kernel_sizes, channels)):
        layers.append(
            torch.nn.Conv1d(
                in_feat, c, kernel_size=k, padding=int((k - 1) / 2)
            )
        )
        layers.append(activation())
        if use_batch_norm:
            layers.append(torch.nn.BatchNorm1d(c))
        in_feat = c

    return torch.nn.Sequential(*layers)


class ReadLevelConv(torch.nn.Module):
    """1-d convs along the read level features."""

    def __init__(
        self,
        in_features=5,
        out_dim=128,
        kernel_sizes=[1, 17],
        channel_dim=128,
        use_batch_norm=True,
    ):
        """Initialise read level conv."""
        super().__init__()
        if isinstance(channel_dim, int):
            channel_dim = [channel_dim] * len(kernel_sizes)

        self.in_features = in_features
        self.out_dim = out_dim
        self.kernel_sizes = kernel_sizes
        self.channel_dim = channel_dim
        self.use_batch_norm = use_batch_norm

        self.convs = make_1dconv_layers(
            kernel_sizes=kernel_sizes,
            in_feat=in_features,
            channels=channel_dim,
            use_batch_norm=use_batch_norm,
        )

        self.expansion_layer = torch.nn.Linear(channel_dim[-1], out_dim)

    def forward(self, x):
        """Forward pass."""
        return self.convs(x)


class MeanPooler(torch.nn.Module):
    """Calculates mean over non-empty read level embeddings."""

    def __init__(self, *args, **kwargs):
        """Initialise mean pooling layer."""
        super().__init__()

    def forward(self, x, non_empty_position_mask):
        """Calculate mean over non-empty read level embeddings.

        :param x: torch.Tensor, read level embeddings,
            shape (batch, read_depth, position_dim, transformer_size)
        :param non_empty_position_mask: torch.Tensor, boolean mask of
            non-empty/padding reads, shape (batch, read_depth)
        """
        read_depths = non_empty_position_mask.sum(-1)
        x = (x * non_empty_position_mask[..., None, None]).sum(
            dim=1
        ) / read_depths[..., None, None]
        return x


class Poolers(enum.Enum):
    """Enum for methods of pooling read embeddings into a single tensor."""

    MEAN = "mean", MeanPooler

    @staticmethod
    def from_str(pooler_type):
        """Map string to Pooler class."""
        for pooler in Poolers:
            if pooler.value[0] == pooler_type:
                return pooler.value[1]
        raise ValueError(f"Unknown PoolerType {pooler_type}")
