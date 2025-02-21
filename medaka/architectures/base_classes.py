"""Base class for all torch medaka models."""
import medaka.features
import medaka.models


class CountsMatrixModel(medaka.models.TorchModel):
    """Base class for models that take counts matrices as input."""

    def get_model_input_features(self, batch):
        """Return the counts matrix from the batch."""
        return batch.counts_matrix

    def check_feature_encoder_compatibility(self, fenc):
        """Check feature encoder is valid for this model."""
        clsname = type(self).__name__

        # check is read level feature encoder
        valid_fencs = [
            medaka.features.CountsFeatureEncoder,
            medaka.features.ReadAlignmentFeatureEncoder
        ]

        if not any(isinstance(fenc, vfenc) for vfenc in valid_fencs):
            msg = f"{type(fenc)} is not a valid feature encoder for {clsname}."
            raise ValueError(msg)


class ReadLevelFeaturesModel(medaka.models.TorchModel):
    """Base class for models that take read level features as input."""

    def get_model_input_features(self, batch):
        """Return the read level features from the batch."""
        features = batch.read_level_features
        return features
