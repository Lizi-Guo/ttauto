# SPDX-License-Identifier: GPL-3.0-or-later
"""Strict, portable experiment specifications (no executable configuration)."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path


def keys(obj, allowed, context):
    if not isinstance(obj, dict):
        raise ValueError(f"{context} must be an object")
    unknown = set(obj) - set(allowed)
    if unknown:
        raise ValueError(f"unknown {context} fields: {sorted(unknown)}")


def positive(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def normalize(data):
    if isinstance(data, dict) and data.get('schema') == 'ttauto-proof-v1':
        from .proof_config import normalize as normalize_proof
        return normalize_proof(data)
    cfg = deepcopy(data)
    keys(cfg, {"schema", "name", "surface", "source", "upper", "pruning", "limits",
               "geometry", "target_polynomial"}, "configuration")
    if cfg.get("schema") != "ttauto-experiment-v1":
        raise ValueError("schema must be ttauto-experiment-v1")
    surface = cfg.get("surface")
    keys(surface, {"kind", "genus", "punctures", "inner_prongs", "outer_prongs",
                   "inner_colors", "foliation"}, "surface")
    if surface.get("kind") not in ("disc", "closed"):
        raise ValueError("surface.kind must be disc or closed")
    genus = surface.setdefault("genus", 0 if surface["kind"] == "disc" else 2)
    positive(genus, "genus", 0)
    if (surface["kind"] == "disc" and genus != 0) or (surface["kind"] == "closed" and genus < 2):
        raise ValueError("discs have genus 0; closed experiments require genus >= 2")
    for field in ("inner_prongs", "outer_prongs"):
        values = surface.get(field)
        if not isinstance(values, list) or not values:
            raise ValueError(f"{field} must be a nonempty list")
        for value in values:
            positive(value, field)
    inner, outer = surface["inner_prongs"], surface["outer_prongs"]
    if sum(p-2 for p in inner+outer) != 4*genus-4:
        raise ValueError("prongs violate the Euler relation")
    colors = surface.setdefault("inner_colors", ["puncture" if surface["kind"] == "disc" and p == 1
                                                else "filled" for p in inner])
    if len(colors) != len(inner) or any(c not in ("puncture", "filled") for c in colors):
        raise ValueError("inner_colors must specify puncture or filled for each inner polygon")
    if any(p < 2 for p, c in zip(inner, colors) if c == "filled"):
        raise ValueError("a one-prong point cannot be filled")
    if surface["kind"] == "disc":
        if len(outer) != 1:
            raise ValueError("a disc has one distinguished outer face")
        n = surface.setdefault("punctures", colors.count("puncture"))
        if positive(n, "punctures", 3) != colors.count("puncture"):
            raise ValueError("puncture count differs from inner_colors")
    elif any(c != "filled" for c in colors) or any(p < 2 for p in outer):
        raise ValueError("closed surfaces fill all faces; one-prong filling is invalid")
    ordered = sorted(zip(inner, colors))
    surface["inner_prongs"] = [p for p, _ in ordered]
    surface["inner_colors"] = [c for _, c in ordered]
    surface["outer_prongs"] = sorted(outer)
    if surface.setdefault("foliation", "any") not in ("any", "orientable", "nonorientable"):
        raise ValueError("unknown foliation type")
    if isinstance(cfg.get("upper"), (float, bool)):
        raise ValueError("upper must be an exact rational string, e.g. '7/4'")
    upper = Fraction(cfg["upper"])
    if upper <= 1:
        raise ValueError("upper must exceed 1")
    cfg["upper"] = str(upper)
    source = cfg.setdefault("source", {"kind": "ribbon"})
    keys(source, {"kind", "path", "executable", "stratum_index"}, "source")
    # The old label is a configuration alias, never a remote provider.
    if source.get("kind") == "upstream":
        source["kind"] = "native"
    if source.get("kind") not in ("ribbon", "automaton", "native"):
        raise ValueError("source.kind must be ribbon, native or automaton")
    if source["kind"] == "native":
        if surface["kind"] != "disc" or not isinstance(source.get("executable"), str):
            raise ValueError("native source requires a disc and a local exporter executable")
        positive(source.setdefault("stratum_index", 0), "stratum_index", 0)
    if source["kind"] == "automaton" and not isinstance(source.get("path"), str):
        raise ValueError("automaton source requires a path")
    pruning = cfg.setdefault("pruning", {})
    keys(pruning, {"norm", "cycle", "homology", "lefschetz", "mod2", "perron", "trace_horizon"}, "pruning")
    for key in ("norm", "cycle", "homology", "lefschetz", "mod2", "perron"):
        if type(pruning.setdefault(key, key == "norm")) is not bool:
            raise ValueError(f"pruning.{key} must be boolean")
    positive(pruning.setdefault("trace_horizon", 12), "trace_horizon")
    limits = cfg.setdefault("limits", {})
    keys(limits, {"events", "depth", "state_decorations"}, "limits")
    for name, value in limits.items():
        positive(value, f"limits.{name}", 0)
    geometry = cfg.setdefault("geometry", {"kind": "builtin"})
    keys(geometry, {"kind"}, "geometry")
    if geometry.get("kind") not in ("audit", "builtin"):
        raise ValueError("geometry.kind must be builtin or audit; external classifiers are not supported")
    polynomial = cfg.get("target_polynomial")
    if polynomial is not None and (not isinstance(polynomial, list) or len(polynomial) < 2
                                  or polynomial[0] != 1 or any(type(c) is not int for c in polynomial)):
        raise ValueError("target_polynomial must be a monic integer polynomial (descending coefficients)")
    return cfg


def load(path):
    path = Path(path).resolve()
    cfg = normalize(json.loads(path.read_text(encoding="utf-8")))
    paths = [(cfg["source"], "path"), (cfg["source"], "executable")]
    if cfg.get('witness'):
        paths.append((cfg['witness'], 'path'))
    for obj, key in paths:
        if key in obj:
            resolved = (path.parent / obj[key]).resolve()
            if key == "executable" and not resolved.exists() and resolved.with_suffix(".exe").exists():
                resolved = resolved.with_suffix(".exe")
            obj[key] = str(resolved)
    return cfg
