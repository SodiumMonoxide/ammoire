"""Tight-binding building blocks and moiré diagnostics for AM/HAM systems."""

from .hyperspin import SquareHAM, charge_channels, projected_block
from .commensurate import COMMENSURATE_PAIRS, CommensurateTriangularBilayer
from .compensated import CompensatedTriangularBilayer
from .models import MODEL_REGISTRY, make_model
from .moire import MoireBilayer

__all__ = [
    "COMMENSURATE_PAIRS", "CommensurateTriangularBilayer",
    "CompensatedTriangularBilayer",
    "MODEL_REGISTRY", "MoireBilayer", "SquareHAM", "charge_channels",
    "make_model", "projected_block",
]
