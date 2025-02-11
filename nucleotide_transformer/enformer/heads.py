import logging
from typing import Dict, List, Optional

import haiku as hk
import jax
import jax.numpy as jnp
from haiku import initializers
from typing_extensions import TypeAlias

from nucleotide_transformer.heads import UNET1DSegmentationHead

SequenceMask: TypeAlias = jnp.ndarray


logger = logging.getLogger(__name__)


class UNetHead(hk.Module):
    """
    Returns a probability between 0 and 1 over a target feature presence
    for each nucleotide in the input sequence. Assumes the sequence has been tokenized
    with non-overlapping k-mers with k being nucl_per_token.
    """

    def __init__(
        self,
        features: List[str],
        num_classes: int = 2,
        embed_dimension: int = 1024,
        nucl_per_token: int = 6,
        num_layers: int = 2,
        remove_cls_token: bool = True,
        name: Optional[str] = None,
    ):
        """
        Args:
            features: List of features names.
            num_classes: Number of classes.
            embed_dimension: Embedding dimension.
            nucl_per_token: Number of nucleotides per token.
            num_layers: Number of layers.
            remove_cls_token: Whether to remove the CLS token.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self._num_features = len(features)
        self._num_classes = num_classes
        self.nucl_per_token = nucl_per_token
        self.remove_cls_token = remove_cls_token

        w_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        b_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        unet = UNET1DSegmentationHead(
            num_classes=embed_dimension // 2,
            output_channels_list=tuple(
                embed_dimension * (2**i) for i in range(num_layers)
            ),
        )
        fc = hk.Linear(
            self.nucl_per_token * self._num_classes * self._num_features,
            w_init=w_init,
            b_init=b_init,
            name="fc",
        )
        self._fc = hk.Sequential([unet, jax.nn.swish, fc])

    def __call__(
        self, x: jnp.ndarray, sequence_mask: SequenceMask
    ) -> Dict[str, jnp.ndarray]:
        """
        Input shape: (batch_size, sequence_length, embed_dim)
        Output_shape: (batch_size, sequence_length * nucl_per_token,
            num_features, num_classes)
        """
        if self.remove_cls_token:
            x = x[:, 1:]
        logits = self._fc(x)
        batch_size, seq_len = x.shape[0], x.shape[1]
        logits = jnp.reshape(
            logits,
            (
                batch_size,
                seq_len * self.nucl_per_token,
                self._num_features,
                self._num_classes,
            ),
        )
        return {"logits": logits}
