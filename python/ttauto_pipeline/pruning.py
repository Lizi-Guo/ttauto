# SPDX-License-Identifier: GPL-3.0-or-later
"""Optional necessary-condition modules. Each has an explicit mathematical scope."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from ttauto_cert.exact_math import charpoly, compare_largest_real_roots, count_real_roots
from ttauto_cert.cycle_pruning import CycleGrowthPruner
from ttauto_cert.lefschetz import compatible_traces, traces
from ttauto_cert.hyperelliptic_mod2 import compatible_branch_cycle_types
from ttauto_cert.perron_obstructions import strict_perron_obstruction
from .topology import orientation_homology


@dataclass(frozen=True)
class Decision:
    outcome: str  # pass, reject, skip; operational exceptions are never rejections
    reason: str
    evidence: object = None


class Context:
    def __init__(self, automaton, config, candidate):
        self.automaton, self.config, self.candidate = automaton, config, candidate
        self.state = automaton.states[candidate["start"]]
        self._words = self._cover = None

    @property
    def mapping(self):
        if self._words is None:
            self._words = self.automaton.compose(self.candidate["start"], self.candidate["path"], self.candidate["terminal"])
        return self._words

    @property
    def cover(self):
        if self._cover is None:
            self._cover = orientation_homology(self.state, *self.mapping)
        return self._cover


class CycleModule:
    name, stage = "cycle", "closed_matrix"
    contract = "completed nonnegative matrix only; never a permuted folding prefix"

    def __init__(self, config):
        self.pruner = CycleGrowthPruner(Fraction(config["upper"]))

    def check(self, context):
        witness = self.pruner.witness(context.candidate["matrix"])
        return Decision("reject" if witness else "pass", "cycle_growth_exceeds_upper" if witness else "no_cycle_obstruction", witness)


class HomologyModule:
    name, stage = "homology", "signed_map"
    contract = "positive tangent-orientation lift on filled H1; formal carrier must realize its prescribed prongs"

    def __init__(self, config):
        pass

    def check(self, context):
        cover = context.cover
        p = cover["polynomial"]
        if len(p) == 1:
            return Decision("reject", "positive_lift_has_no_expanding_homology", cover)
        if p != p[::-1]:
            raise ValueError("orientation-cover homology is not reciprocal: invalid map or implementation")
        if not count_real_roots(p) or compare_largest_real_roots(p, context.candidate["polynomial"]) != 0:
            return Decision("reject", "positive_lift_growth_differs_from_real_transition", cover)
        return Decision("pass", "positive_lift_growth_matches", cover)


class LefschetzModule:
    name, stage = "lefschetz", "signed_map"
    contract = "finite-horizon necessary periodic-orbit conditions; both lifts tested on a connected orientation cover"

    def __init__(self, config):
        self.horizon = config["pruning"]["trace_horizon"]

    def check(self, context):
        cover = context.cover
        if cover["genus"] < 1:
            return Decision("skip", "cover_has_no_positive_genus")
        groups = sorted(Counter(p for p in cover["lifted_prongs"] if p >= 4).items())
        if cover["connected_double_cover"]:
            lifts = [(cover["polynomial"], False), (cover["other_lift_polynomial"], True)]
        else:
            lifts = [(charpoly(cover["actual_base_matrix"]), cover["base_action_sign"] < 0)]
        signs = [negative for _, negative in lifts]
        counts = [len(compatible_traces(traces(p, self.horizon)[1:], len(p)-1,
                                       groups, negative=negative)) for p, negative in lifts]
        evidence = {"horizon": self.horizon, "negative_lifts": list(signs), "compatible_patterns": counts,
                    "groups": groups, "actual_lift_polynomials": [p for p, _ in lifts],
                    "deck_matrix": cover.get("deck_matrix")}
        return Decision("reject" if 0 in counts else "pass", "lefschetz_orbit_conditions", evidence)


class Mod2Module:
    name, stage = "mod2", "signed_map"
    contract = "sphere double cover only; colored odd-prong branch permutations"

    def __init__(self, config):
        pass

    def check(self, context):
        surface = context.config["surface"]
        if surface["kind"] != "disc":
            return Decision("skip", "base_is_not_a_sphere")
        cover = context.cover
        if not cover["connected_double_cover"]:
            return Decision("skip", "orientation_cover_is_disconnected")
        groups = list(Counter((p, c) for p, c in zip(surface["inner_prongs"], surface["inner_colors"]) if p % 2).values())
        fixed = surface["outer_prongs"][0] % 2
        matches = compatible_branch_cycle_types(cover["polynomial"], groups, fixed)
        return Decision("pass" if matches else "reject", "colored_branch_permutation_mod2",
                        {"groups": groups, "fixed": fixed, "compatible_types": len(matches)})


class PerronModule:
    name, stage = "perron", "signed_map"
    contract = "simple strictly dominant positive eigenvalue on the orienting cover"

    def __init__(self, config):
        pass

    def check(self, context):
        p = context.cover["polynomial"]
        from ttauto_cert.exact_math import compare_largest_real_root
        if len(p) == 1 or not count_real_roots(p) or compare_largest_real_root(p, 1) <= 0:
            return Decision("reject", "no_positive_expanding_lift_root")
        result = strict_perron_obstruction(p)
        return Decision("reject" if result["excluded"] else "pass", result["reason"], result)


REGISTRY = {cls.name: cls for cls in (CycleModule, HomologyModule, Mod2Module, LefschetzModule, PerronModule)}


def modules(config):
    return [cls(config) for name, cls in REGISTRY.items() if config["pruning"][name]]
