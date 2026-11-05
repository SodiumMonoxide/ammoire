"""Optional adapters to the models that were validated in the earlier reports."""

from dataclasses import dataclass
from importlib import import_module
import os
from pathlib import Path
import sys

import numpy as np


# module, function, basis dimension, length in the source convention, unit
SOURCES = {
    "roig-2d": ("legacy_bridge", "minimal2d.H", 4, 1.0, "source model units"),
    "farajollahpour": ("legacy_bridge", "c4k2d.H", 2, 1.0, "eV"),
    "orbital-hall-6": ("legacy_bridge", "legacy6.H", 6, 1.0, "eV"),
    "wen-kee": ("meta_altermagnet", "H_primitive", 4, 1.0, "|t|"),
    "gese-90": ("moire_gese", "H", 2, 25.658, "meV"),
    "guo-90": ("twisted_bilayer", "fourier_90", 4, 1.0, "t"),
}


@dataclass(frozen=True)
class _DailyModel:
    name: str
    size: int
    a: float
    energy_unit: str
    matrix: object
    parameters: dict
    energy_scale: float
    dimension = 2
    geometry = "square"

    @property
    def lattice(self):
        return self.a * np.eye(2)

    @property
    def reciprocal(self):
        return 2 * np.pi / self.a * np.eye(2)

    def hamiltonian(self, k):
        k = np.asarray(k, dtype=float)
        if k.ndim == 0 or k.shape[-1] != 2:
            raise ValueError("k must have shape (..., 2)")
        x, y = np.moveaxis(k * self.a, -1, 0)
        h = self.matrix(x, y, **self.parameters)
        if self.name in {"roig-2d", "wen-kee"}:
            # Cell gauge for an A/B embedding at (0,0), (a/2,a/2).
            phase = np.exp(0.5j * (x + y))
            u = np.stack((np.ones_like(phase), np.ones_like(phase), phase, phase), axis=-1)
            h = u[..., :, None].conj() * h * u[..., None, :]
        return self.energy_scale * h


def daily_model(name, *, continuation=None, a=None, energy_scale=1.0,
                energy_unit=None, **parameters):
    """Reuse the daily-report code through its original bridge.

    Set ALTERMAGNET_TB_ROOT to the older project if it has moved. Guo-90
    and GeSe-90 are completed bilayer benchmarks to use as standalone bands.
    """
    module, symbol, size, default_a, unit = SOURCES[name]
    if continuation is None:
        continuation = os.environ.get("ALTERMAGNET_CONTINUATION")
    if continuation is None:
        raise FileNotFoundError(
            "Set ALTERMAGNET_CONTINUATION or pass continuation= for the local daily model library"
        )
    continuation = Path(continuation).resolve()
    if not (continuation / "tb_extensions").is_dir():
        raise FileNotFoundError(f"tb_extensions missing in {continuation}")
    sys.path.insert(0, str(continuation))
    try:
        matrix = import_module(f"tb_extensions.{module}")
    finally:
        sys.path.pop(0)
    for part in symbol.split("."):
        matrix = getattr(matrix, part)
    a = default_a if a is None else float(a)
    if not np.isfinite(a) or a <= 0 or not np.isfinite(energy_scale) or energy_scale <= 0:
        raise ValueError("a and energy_scale must be positive and finite")
    if name == "gese-90":
        field = parameters.pop("field", 0.0)
        if parameters or field not in (0.0, 0.1):
            raise ValueError("GeSe preset field must be 0.0 or 0.1 eV/Å")
        source = import_module("tb_extensions.moire_gese")
        hoppings = source.FIELD_01 if field == 0.1 else source.ZERO_FIELD
        matrix = lambda x, y: source.H(x, y, *hoppings)
    return _DailyModel(name, size, a, energy_unit or unit, matrix, parameters, energy_scale)
