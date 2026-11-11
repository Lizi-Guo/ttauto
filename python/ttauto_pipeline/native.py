# SPDX-License-Identifier: GPL-3.0-or-later
"""Build automata with the in-tree C++ signed-fold algorithm."""
import json
from pathlib import Path
import subprocess
import tempfile
from .combinatorics import require
from .ribbon import real_matrix
from .ribbon_states import graph_from_fans


def convert(data, surface):
    require(data.get("schema") in ("ttauto-native-export-v1", "ttauto-upstream-export-v1"), "native export schema")
    require(data["punctures"] == surface["punctures"], "native punctures")
    states, renamings, orders = [], [], []
    for row in data["states"]:
        m, remaining, polygons = row["main_edges"], set(range(len(row["fans"]))), []
        while remaining:
            start = min(remaining)
            at, vertices = start, []
            while at in remaining:
                vertices.append(at)
                remaining.remove(at)
                at = row["side_to"][at]
            require(at == start, "native side permutation")
            colors = {row["punctured"][v] for v in vertices}
            require(len(colors) == 1, "native polygon puncturing")
            polygons.append((len(vertices), "puncture" if colors.pop() else "filled", vertices))
        polygons.sort()
        order = [v for _, _, vertices in polygons for v in vertices]
        naming = {e: e for e in range(1, m+1)}
        naming.update({m+old+1: m+new+1 for new, old in enumerate(order)})
        graph = graph_from_fans([p for p, _, _ in polygons], [row["fans"][v] for v in order])
        require(graph is not None, "disconnected native graph")
        graph["inner_colors"] = [c for _, c, _ in polygons]
        states.append({"graph": graph, "arcs": []})
        renamings.append(naming)
        orders.append(order)
    for si, row in enumerate(data["states"]):
        source_naming, source_order = renamings[si], orders[si]
        for arc in row["arcs"]:
            target = arc["target"]
            naming, order = renamings[target], orders[target]
            words = [None]*len(source_naming)
            for e, word in enumerate(arc["edge_images"], 1):
                words[source_naming[e]-1] = [naming[abs(d)]*(1 if d > 0 else -1) for d in word]
            vertices = [order.index(arc["vertex_images"][v]) for v in source_order]
            states[si]["arcs"].append({"target": target, "vertex_images": vertices,
                                      "edge_images": words, "matrix": real_matrix(words, row["main_edges"])})
    return {"schema": "ttauto-automaton-v1", "surface": surface, "states": states,
            "state_inventory_complete": False,
            "coverage": "native fold component from selected seed; full state coverage is not inferred from reachability",
            "enumeration": {"provider": "in-tree-cpp", "stratum_index": data["stratum_index"]}}


def construct(config):
    source = config["source"]
    with tempfile.TemporaryDirectory(prefix="ttauto-export-") as directory:
        path = Path(directory)/"graph.json"
        subprocess.run([source["executable"], str(config["surface"]["punctures"]),
                        str(source["stratum_index"]), str(path)], check=True,
                       capture_output=True, text=True,
                       creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        return convert(json.loads(path.read_text(encoding="utf-8")), config["surface"])
