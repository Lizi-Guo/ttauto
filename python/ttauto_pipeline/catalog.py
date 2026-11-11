# SPDX-License-Identifier: GPL-3.0-or-later
"""Enumerate experiment specifications, never insert tabulated minima."""
from .config import normalize
from ttauto_cert.lefschetz import partitions


def disc_strata(n, upper, one_prong=False):
    if type(n) is not int or n < 3:
        raise ValueError("n must be >= 3")
    for outer in range(1, n-1):
        budget = n-2-outer
        for puncture_budget in ([0] if one_prong else range(budget+1)):
            for increments in partitions(puncture_budget):
                if len(increments) > n:
                    continue
                punctures = [1]*(n-len(increments))+[1+x for x in increments]
                for zero_orders in partitions(budget-puncture_budget):
                    zeros = [2+x for x in zero_orders]
                    yield normalize({"schema": "ttauto-experiment-v1",
                                     "name": f"D{n}-P{'-'.join(map(str,punctures))}-Z{'-'.join(map(str,zeros))}-O{outer}",
                                     "surface": {"kind": "disc", "punctures": n,
                                                 "inner_prongs": punctures+zeros,
                                                 "inner_colors": ["puncture"]*n+["filled"]*len(zeros),
                                                 "outer_prongs": [outer]}, "upper": upper})


def closed_strata(genus, upper, auxiliary_points=1):
    if type(genus) is not int or genus < 2 or type(auxiliary_points) is not int or auxiliary_points < 1:
        raise ValueError("genus >= 2 and auxiliary_points >= 1 required")
    for kind, degree in (("H", 2*genus-2), ("Q", 4*genus-4)):
        for orders in partitions(degree):
            prongs = [2*d+2 if kind == "H" else d+2 for d in orders]
            yield normalize({"schema": "ttauto-experiment-v1", "name": f"g{genus}-{kind}{'-'.join(map(str,orders))}-regular{auxiliary_points}",
                             "surface": {"kind": "closed", "genus": genus,
                                         "inner_prongs": [2]*auxiliary_points, "outer_prongs": prongs,
                                         "foliation": "orientable" if kind == "H" else "nonorientable"}, "upper": upper})
