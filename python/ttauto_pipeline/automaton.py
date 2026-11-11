# SPDX-License-Identifier: GPL-3.0-or-later
"""A single finite folding graph interface for native and Python state builders."""
from copy import deepcopy
import json
from pathlib import Path
from ttauto_cert.exact_math import identity, matrix_multiply
from ttauto_cert.standard_maps import exact_sum_bound
from ttauto_cert.ribbon_verify import substitute, inverse_elementary, cyclic_word
from .combinatorics import require
from .ribbon import SurfaceState, build, real_matrix
from .topology import terminals, verify_map, ribbon_graph
from .ribbon_states import graph_from_fans
from .carrying import certify
from .carrying_verify import verify as verify_carrying


class Automaton:
    def __init__(self, data, surface):
        require(data.get("schema") == "ttauto-automaton-v1", "automaton schema")
        require(data.get("surface") == surface, "automaton surface differs from experiment")
        self.data = data
        self.states = [SurfaceState(row["graph"]) for row in data["states"]]
        require(bool(self.states), "empty automaton")
        self.arcs = [row["arcs"] for row in data["states"]]
        self.terminals = []
        self.dimension = self.states[0].m
        for state in self.states:
            graph = state.graph
            require(type(state.m) is int and state.m > 0, "invalid real dimension")
            require(all(type(d) is int and d != 0 for fan in graph["real_fans"] for d in fan), "nonintegral real darts")
            rebuilt = graph_from_fans(graph["inner_prongs"], graph["real_fans"])
            require(rebuilt is not None and all(rebuilt[key] == graph[key] for key in
                    ("main_edges", "genus", "cyclic_directions", "boundary_cycles", "outer_prongs")), "inconsistent ribbon state")
        for index, (state, row) in enumerate(zip(self.states, self.arcs)):
            require(state.m == self.dimension, "variable real dimension")
            require(state.graph["genus"] == surface["genus"] and
                    state.graph["inner_prongs"] == surface["inner_prongs"] and
                    state.graph["outer_prongs"] == surface["outer_prongs"] and
                    state.graph["inner_colors"] == surface["inner_colors"], "state topological type")
            for arc in row:
                require(type(arc["target"]) is int and 0 <= arc["target"] < len(self.states), "arc target")
                target = self.states[arc["target"]]
                require(all(type(v) is int for v in arc["vertex_images"]), "nonintegral vertex map")
                verify_map(state, target, arc["vertex_images"], arc["edge_images"])
                source_polygons, target_polygons = ribbon_graph(state.graph), ribbon_graph(target.graph)
                targets = {cyclic_word(w): i for i, w in enumerate(target_polygons["polygons"])}
                for pi, polygon in enumerate(source_polygons["polygons"]):
                    key = cyclic_word(substitute(arc["edge_images"], polygon))
                    require(key in targets, "inner polygon maps to outer face")
                    qi = targets[key]
                    require((state.graph["inner_prongs"][pi], state.graph["inner_colors"][pi]) ==
                            (target.graph["inner_prongs"][qi], target.graph["inner_colors"][qi]), "polygon color transport")
                inverse_elementary(arc["edge_images"])
                matrix = real_matrix(arc["edge_images"], state.m)
                require(all(type(x) is int for row in arc["matrix"] for x in row) and matrix == arc["matrix"], "signed word/matrix mismatch")
                require(sum(map(sum, matrix)) == state.m+1, "not an elementary expanding fold")
                # Every column and row must contain a permutation support. A
                # positive row sum alone would not prove monotone norm growth.
                one = [sum(row) for row in matrix]
                require(sorted(one) == [1]*(state.m-1)+[2], "fold row counts")
                columns = [sum(row[j] for row in matrix) for j in range(state.m)]
                require(sorted(columns) == [1]*(state.m-1)+[2] and
                        matrix[one.index(2)][columns.index(2)] > 0, "fold lacks permutation-plus-one support")
                cert = certify(state.graph, target.graph, arc)
                orders = cert["positive_band_lane_orders"]
                verify_carrying(state.graph, target.graph, arc,
                                [orders[e] for e in range(1, state.e+1)])
            self.terminals.append(terminals(state))

    @classmethod
    def construct(cls, config):
        surface = config["surface"]
        if config["source"]["kind"] == "automaton":
            data = json.loads(Path(config["source"]["path"]).read_text(encoding="utf-8"))
        elif config["source"]["kind"] in ("native", "upstream"):
            from .native import construct
            data = construct(config)
        else:
            states, arcs, counts = build(surface["inner_prongs"], surface["outer_prongs"],
                                        surface["genus"], surface["inner_colors"],
                                        config["limits"].get("state_decorations"))
            data = {"schema": "ttauto-automaton-v1", "surface": surface,
                    "enumeration": counts, "state_inventory_complete": True,
                    "coverage": "prescribed inner/outer type; not all types on a closed surface",
                    "states": [{"graph": s.graph, "arcs": row} for s, row in zip(states, arcs)]}
        return cls(data, surface)

    def compose(self, start, path, terminal):
        require(type(start) is int and 0 <= start < len(self.states), "start state")
        require(type(terminal) is int and 0 <= terminal < len(self.terminals[start]), "terminal index")
        source = self.states[start]
        words = [[e] for e in range(1, source.e+1)]
        vertices = list(range(source.k))
        at = start
        for index in path:
            require(type(index) is int and 0 <= index < len(self.arcs[at]), "path arc index")
            arc = self.arcs[at][index]
            words = [substitute(arc["edge_images"], w) for w in words]
            vertices = [arc["vertex_images"][v] for v in vertices]
            at = arc["target"]
        require(at == start and path, "nonempty closed path required")
        closure = self.terminals[start][terminal]
        words = [substitute(closure["edge_images"], w) for w in words]
        vertices = [closure["vertex_images"][v] for v in vertices]
        verify_map(source, source, vertices, words)
        return words, vertices


class FoldSource:
    """Resumable DFS without matrix-only deduplication or powers of maps.

    Cursor contains the actual frontier. Only elementary permutation-plus-one
    folds are admitted. Total entry sum increases by at least one per fold;
    the Ham--Song bound therefore also gives a finite depth bound. Terminal
    symmetries are ALL oriented, color-preserving ribbon automorphisms.
    """
    def __init__(self, automaton, config, cursor=None):
        self.automaton = automaton
        self.bound = exact_sum_bound(config["upper"], automaton.dimension)
        self.theoretical_depth = self.bound-automaton.dimension
        self.depth = min(config["limits"].get("depth", self.theoretical_depth), self.theoretical_depth)
        self.depth_truncated = self.depth < self.theoretical_depth
        self.cursor = deepcopy(cursor) if cursor else {"next_start": 0, "stack": []}

    def reject_prefix(self):
        require(bool(self.cursor["stack"]), "no active prefix")
        self.cursor["stack"].pop()

    def next_event(self):
        a, c = self.automaton, self.cursor
        while True:
            if not c["stack"]:
                if c["next_start"] >= len(a.states):
                    return None
                start = c["next_start"]
                c["next_start"] += 1
                c["stack"].append({"start": start, "state": start, "path": [],
                                    "matrix": identity(a.dimension), "arc": 0, "terminal": 0})
            node = c["stack"][-1]
            if node["path"] and node["state"] == node["start"] and node["terminal"] < len(a.terminals[node["start"]]):
                ti = node["terminal"]
                node["terminal"] += 1
                return {"kind": "candidate", "start": node["start"], "path": node["path"][:],
                        "terminal": ti, "matrix": matrix_multiply(a.terminals[node["start"]][ti]["matrix"], node["matrix"])}
            if len(node["path"]) >= self.depth or node["arc"] == len(a.arcs[node["state"]]):
                c["stack"].pop()
                continue
            ai = node["arc"]
            node["arc"] += 1
            arc = a.arcs[node["state"]][ai]
            matrix = matrix_multiply(arc["matrix"], node["matrix"])
            require(sum(map(sum, matrix)) > sum(map(sum, node["matrix"])), "nonexpanding prefix")
            child = {"start": node["start"], "state": arc["target"], "path": node["path"]+[ai],
                     "matrix": matrix, "arc": 0, "terminal": 0}
            c["stack"].append(child)
            return {"kind": "prefix", "matrix": matrix, "depth": len(child["path"])}
