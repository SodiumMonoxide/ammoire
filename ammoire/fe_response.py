"""Relate a spin Laue prediction to a spin-resolved two-domain band test."""

import json
from pathlib import Path

import numpy as np

from .fe_reference import literature_cases
from .fe_spectra import audit_domain_spectra, audit_group_spectra
from .fe_symmetry import SpinLaueGroup
from .models import MomentumSlice, make_model


def domain_response(group, operation, positive, negative, kpoints, *, spin_axis=(0, 0, 1),
                    atol=1e-8):
    """Test a specified switching path between two collinear TB Hamiltonians.

    The caller supplies the full signed Laue group of the positive domain
    and the magnetic-configuration-preserving map from positive to negative.
    The returned class is conditional on these inputs. Residuals measure
    absolute energy differences on the supplied Cartesian momenta.
    """
    if group.dimension != positive.dimension or group.dimension != negative.dimension:
        raise ValueError("group and both models must have the same momentum dimension")
    if any(isinstance(model, MomentumSlice) and model.kz != 0 for model in (positive, negative)):
        raise ValueError("use the 3D bulk models for sections at nonzero kz")
    response = group.classify(operation)
    action = group.momentum_action(operation)
    domain_test = audit_domain_spectra(positive, negative, kpoints, action,
                                       spin_axis=spin_axis, atol=atol)
    group_test = audit_group_spectra(positive, kpoints, group, spin_axis=spin_axis, atol=atol)
    relation = domain_test["same_k_relation"]
    required_residual = {"I": "same_k_preserved_residual", "II": "same_k_reversed_residual"}
    if not group_test["matches"] or not domain_test["mapped_same_spin_matches"]:
        verdict = "inconsistent"
    elif response in required_residual and domain_test[required_residual[response]] > domain_test["atol"]:
        verdict = "inconsistent"
    elif relation == "unresolved" or (response == "III" and relation != "neither"):
        verdict = "undistinguished_on_sample"
    else:
        verdict = "consistent_on_sample"
    return {"response_type": response, "sample_verdict": verdict,
            "positive_domain_group": group_test, "two_domain_spectra": domain_test}


def _registry_model(spec):
    return make_model(spec["model"], **spec.get("parameters", {}))


def bilayer_demo():
    """Build three switched structures and compare their spin spectra."""
    positive = make_model("sun-hex-bilayer")
    mirror_z = np.diag([1, 1, -1])
    mirror_y = np.diag([1, -1, 1])
    sine = np.sqrt(3) / 2
    c3z = np.array([[-0.5, -sine, 0], [sine, -0.5, 0], [0, 0, 1]])
    group = SpinLaueGroup(2, [(mirror_y, -1)])
    # Both generic and nodal points, expressed in reciprocal primitive coordinates.
    fractions = [[0, 0], [0.17, 0.31], [0.28, -0.19], [-0.43, 0.07], [0.5, 0]]
    points = fractions @ positive.reciprocal
    reports = []
    for label, op, stacking, sign in (
        ("Mz", mirror_z, 3, -1),
        ("Mz My", mirror_z @ mirror_y, 3, 1),
        ("Mz C3z", mirror_z @ c3z, 5, -1),
    ):
        negative = make_model("sun-hex-bilayer", stacking=stacking, layer_sign=sign)
        reports.append({"operation": label, "stacking": stacking, "layer_sign": sign,
                        **domain_response(group, op, positive, negative, points)})
    return {"model": positive.name, "kpoints_cartesian": points.tolist(), "paths": reports}


def _read_path(path):
    spec = json.loads(path.read_text(encoding="utf-8-sig"))
    group = SpinLaueGroup(spec["dimension"],
                          [(g["matrix"], g["spin_sign"]) for g in spec["generators"]])
    op = spec["operation"]
    if "positive" not in spec and "negative" not in spec:
        if "kpoints" in spec:
            raise ValueError("kpoints require both positive and negative domain models")
        return {"response_type": group.classify(op), "group_order": len(group.operations),
                "momentum_action": group.momentum_action(op).tolist()}
    return domain_response(group, op, _registry_model(spec["positive"]),
                           _registry_model(spec["negative"]), spec["kpoints"],
                           spin_axis=spec.get("spin_axis", (0, 0, 1)), atol=spec.get("atol", 1e-8))


def add_command(commands):
    parser = commands.add_parser("fe-response", help="classify and test ferroelectric switching paths")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("manifest", type=Path, nargs="?", help="JSON symmetry and domain specification")
    source.add_argument("--bilayer", action="store_true", help="run the three eight-band TB domain tests")
    parser.add_argument("--out", type=Path, help="also write the JSON report to a file")
    return parser


def run(args, parser):
    try:
        if args.manifest is not None:
            report = _read_path(args.manifest)
        elif args.bilayer:
            report = bilayer_demo()
        else:
            report = literature_cases()
        output = json.dumps(report, ensure_ascii=True, allow_nan=False, indent=2)
        if args.out is not None:
            if args.manifest is not None and args.out.resolve() == args.manifest.resolve():
                raise ValueError("choose an output distinct from the input manifest")
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(output + "\n", encoding="utf-8")
    except (KeyError, ValueError, OSError, TypeError) as error:
        parser.error(str(error))
    print(output)
