"""Read pileup model that takes the majority vote of the pileup.

This model exists to provide a simple baseline for comparison with more
sophisticated models. It is not expected to perform well on real data, nor
to be used outside of testing and debugging.
"""

import warnings

import torch

import medaka.architectures.base_classes as base_classes
import medaka.common


class MajorityVoteModel(base_classes.CountsMatrixModel):
    """Implementation of bidirectional GRU model in pytorch."""

    def __init__(self, time_steps=None, **kwargs):
        """Initialise simple consensus model taking majority vote of pileup."""
        super().__init__()
        self.num_classes = 5  # needed for tests to pass

        if time_steps is not None:
            warnings.warn("timesteps is no longer required to be specified")

        # this is a dummy parameter to ensure that .backward() can be called
        # on this model without torch throwing an error.
        self.dummy_parameter = torch.nn.Parameter(
            torch.zeros(1, requires_grad=True)
        )
        self.config = {
            "model_type": "majority_vote",
            "model_args": {**kwargs},
        }

    def forward(self, pileup, **kwargs):
        """Model forward pass.

        :param x: torch.Tensor, input tensor, shape (batch, n_pos, features).
        :param candidate_positions: batch-length list of lists of candidate
            positions for each read.

        :returns: torch.Tensor, shape (num_candidatees, classes).
        """
        b2i = medaka.common.base2index
        b = pileup[..., b2i['a']:b2i['t']+1] + pileup[..., b2i['A']:b2i['T']+1]
        # sum deletion counts (indexing in this way retains correct shape)
        d = pileup[..., b2i['d']:b2i['d']+1] + pileup[..., b2i['D']:b2i['D']+1]
        # d first as it is the deletion class, labelled 0 in training data
        out = torch.cat([d, b], -1)
        out[..., 0] += 1 - out.sum(-1)
        return out + 0 * self.dummy_parameter
