"""Input metadata and Wigner–Seitz shifts for Wannier90 interpolation."""

from pathlib import Path
import re

import numpy as np


_BOHR_TO_ANGSTROM = 0.529177210903


def seed_file(seed, suffix):
    seed = Path(seed)
    return seed.with_name(seed.name + suffix)


def read_win(path):
    """Read the row-vector cell and interpolation settings from seed.win."""
    lines = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        clean = re.split(r"[!#]", line, maxsplit=1)[0].strip().lower()
        if clean:
            lines.append(clean)
    starts = [i for i, line in enumerate(lines) if line.split() == ["begin", "unit_cell_cart"]]
    ends = [i for i, line in enumerate(lines) if line.split() == ["end", "unit_cell_cart"]]
    if len(starts) != 1 or len(ends) != 1 or ends[0] <= starts[0]:
        raise ValueError(f"{path.name} must contain one unit_cell_cart block")
    block = lines[starts[0] + 1:ends[0]]
    units = {"bohr", "ang", "angstrom", "angstroms"}
    unit = block.pop(0) if block and block[0] in units else "ang"
    try:
        lattice = np.array([[float(x.replace("d", "e")) for x in row.replace(",", " ").split()]
                            for row in block])
    except ValueError as exc:
        raise ValueError(f"invalid lattice in {path.name}") from exc
    if lattice.shape != (3, 3):
        raise ValueError(f"{path.name} unit_cell_cart must have three Cartesian rows")
    if unit == "bohr":
        lattice *= _BOHR_TO_ANGSTROM
    settings = {}
    for line in lines:
        fields = re.split(r"[\s=:]+", line)
        name = fields[0]
        if name not in {"use_ws_distance", "write_ndegen_applied"}:
            continue
        if name in settings or len(fields) != 2 or fields[1].strip(".") not in {"true", "false", "t", "f"}:
            raise ValueError(f"invalid or repeated {name} in {path.name}")
        settings[name] = fields[1].strip(".") in {"true", "t"}
    return lattice, settings


def _header_flag(header, name):
    match = re.search(rf"{name}=\.(true|false)\.", header)
    return match[1] == "true" if match else None


def interpolation_settings(path, settings):
    """Resolve the file header and .win flags before applying shifts."""
    sidecar = None
    header_ws = header_expanded = None
    if path.name.endswith("_hr.dat"):
        candidate = path.with_name(path.name.removesuffix("_hr.dat") + "_wsvec.dat")
        if candidate.is_file():
            sidecar = candidate
            with sidecar.open(encoding="utf-8-sig") as file:
                header = "".join(file.readline().lower().split())
            header_ws = _header_flag(header, "use_ws_distance")
            header_expanded = _header_flag(header, "write_ndegen_applied")
            if header_ws is None and header_expanded is not True:
                raise ValueError(f"{sidecar.name} has no recognized Wigner–Seitz settings")
    expected_ws = settings.get("use_ws_distance") if settings is not None else None
    expected_expanded = settings.get("write_ndegen_applied") if settings is not None else None
    for flag, expected, observed in (("use_ws_distance", expected_ws, header_ws),
                                      ("write_ndegen_applied", expected_expanded, header_expanded)):
        if expected is not None and observed is not None and expected != observed:
            raise ValueError(f"{flag} differs between .win and {sidecar.name}")
    if header_expanded is True or expected_expanded is True:
        if sidecar is not None and header_expanded is not True:
            raise ValueError(f"{sidecar.name} does not mark the expanded export specified by .win")
        return "expanded", None
    use_ws = expected_ws if expected_ws is not None else header_ws
    if use_ws is None:
        # Seed loading uses the Wannier90 3.1 default. The legacy hr-only
        # loader relies on the caller knowing the export mode without a sidecar.
        use_ws = settings is not None
    if use_ws and sidecar is None:
        raise ValueError("Wigner–Seitz interpolation requires the matching *_wsvec.dat sidecar")
    return ("wsvec", sidecar) if use_ws else ("degeneracy", None)


def _ws_integers(line, size):
    fields = line.split() if line is not None else []
    if line is not None and len(fields) != size and len(line) == 5 * size:
        fields = [line[i:i + 5] for i in range(0, len(line), 5)]
    try:
        if len(fields) != size:
            raise ValueError
        return tuple(map(int, fields))
    except ValueError as exc:
        raise ValueError("invalid or truncated Wigner–Seitz shift record") from exc


def expand_wsvec(path, translations, hopping):
    """Accumulate already R-weighted hoppings at R+T, averaging over T.

    Source: https://github.com/wannier-developers/wannier90/blob/v3.1.0/src/plot.F90
    """
    size = hopping.shape[1]
    indices = {r: i for i, r in enumerate(translations)}
    seen = np.zeros(hopping.shape, dtype=bool)
    effective = {}
    with path.open(encoding="utf-8-sig") as file:
        file.readline()
        rows = (line.rstrip("\r\n") for line in file if line.strip())
        for _ in range(hopping.size):
            x, y, z, m, n = _ws_integers(next(rows, None), 5)
            i = indices.get((x, y, z))
            if i is None or not (1 <= m <= size and 1 <= n <= size):
                raise ValueError("Wigner–Seitz shift does not match the hr translations or orbital count")
            m, n = m - 1, n - 1
            if seen[i, m, n]:
                raise ValueError("duplicate Wigner–Seitz shift entry")
            seen[i, m, n] = True
            count, = _ws_integers(next(rows, None), 1)
            if count < 1:
                raise ValueError("Wigner–Seitz shift multiplicity must be positive")
            shifts = set()
            weighted = hopping[i, m, n] / count
            for _ in range(count):
                shift = _ws_integers(next(rows, None), 3)
                if shift in shifts:
                    raise ValueError("duplicate vector in Wigner–Seitz shift list")
                shifts.add(shift)
                r = (x + shift[0], y + shift[1], z + shift[2])
                if r not in effective:
                    effective[r] = np.zeros((size, size), dtype=complex)
                effective[r][m, n] += weighted
        if next(rows, None) is not None:
            raise ValueError("extra records after the expected Wigner–Seitz shift blocks")
    vectors = sorted(effective)
    return vectors, np.array([effective[r] for r in vectors])
