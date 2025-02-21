"""Bidirectional GRU on counts matrix."""

import warnings

import torch

import medaka.architectures.base_classes as base_classes


class GRUModel(base_classes.CountsMatrixModel):
    """Implementation of bidirectional GRU model in pytorch."""

    def __init__(
        self,
        num_features=10,
        num_classes=5,
        gru_size=128,
        n_layers=2,
        bidirectional=True,
        time_steps=None,
        classify_activation=None,
    ):
        """Initialise bidirectional gru model.

        :param num_features: int, number of features for each pileup column.
        :param num_classes: int, number of output class labels.
        :param gru_size: int, size of each gru layer (in each direction).
        :param n_layers: int, number of gru layers.
        :param classify_activation: str, activation to use in classification
        :param timesteps: dummy argument, no longer used
        """
        super().__init__()

        if time_steps is not None:
            warnings.warn("timesteps is no lnoger required to be specified")

        if classify_activation is not None:
            warnings.warn("classify_activation is no longer used")

        self.gru_size = gru_size
        self.num_classes = num_classes
        self.num_features = num_features
        self.n_layers = n_layers
        self.bidirectional = bidirectional

        self.gru = torch.nn.GRU(
            num_features,
            gru_size,
            num_layers=n_layers,
            bidirectional=bidirectional,
            batch_first=True,
        )
        self.linear = torch.nn.Linear(
            2 * gru_size if bidirectional else gru_size, 5
        )
        self.normalise = True

    def forward(self, x):
        """Model forward pass.

        :param x: torch.Tensor, counts matrix input tensor, shape
            (batch size, positions, 2*alphabet size)
        :return: torch.Tensor, predictions for each position,
            shape (batch size, positions, num_classes)
        """
        x = self.gru(x)[0]
        x = self.linear(x)
        if self.normalise:
            # switch normalise off for training to return logits for
            # cross-entropy loss
            x = torch.softmax(x, dim=-1)
        return x
