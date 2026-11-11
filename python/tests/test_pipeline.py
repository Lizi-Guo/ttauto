# SPDX-License-Identifier: GPL-3.0-or-later
from copy import deepcopy
from contextlib import closing
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from ttauto_pipeline.automaton import Automaton, FoldSource
from ttauto_pipeline.catalog import disc_strata, closed_strata
from ttauto_pipeline.config import normalize
from ttauto_pipeline.engine import evaluate, run
from ttauto_pipeline.pruning import Context, modules
from ttauto_pipeline.ribbon import real_matrix
from ttauto_pipeline.ribbon_states import StateBudgetExceeded
from ttauto_pipeline.storage import digest
from ttauto_pipeline.verify import verify

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).with_name("fixtures")


def disc_config(**extra):
    config = dict(schema="ttauto-experiment-v1", surface=dict(kind="disc", inner_prongs=[1, 1, 1], outer_prongs=[1]), upper="8/3")
    config.update(extra)
    return normalize(config)


def records(folder):
    with closing(sqlite3.connect(Path(folder)/"run.sqlite3")) as db:
        return [json.loads(row[0]) for row in db.execute("SELECT payload FROM candidates ORDER BY id")]


def survivors(rows):
    return {(r["start"], tuple(r["path"]), r["terminal"], tuple(r["polynomial"]))
            for r in rows if r["status"] in ("unresolved", "attained")}


class ConfigTests(unittest.TestCase):
    def test_legacy_provider_alias_uses_local_native_algorithm(self):
        native = disc_config(source=dict(kind="native", executable="build/bin/ttauto_pipeline_export"))
        legacy = disc_config(source=dict(kind="upstream", executable="build/bin/ttauto_pipeline_export"))
        self.assertEqual(native, legacy)
        self.assertEqual(native["source"]["kind"], "native")

    def test_exact_input_and_typo_rejection(self):
        for changes in ({"upper": 1.75}, {"pruning": {"homlogy": True}}, {"limits": {"events": -1}}):
            with self.assertRaises(ValueError):
                disc_config(**changes)
        with self.assertRaises(ValueError):
            disc_config(surface=dict(kind="closed", genus=2, inner_prongs=[2], outer_prongs=[8]))

    def test_catalog_euler_and_known_paper_counts(self):
        self.assertEqual(len(list(disc_strata(6, "2", one_prong=True))), 7)
        self.assertEqual(len(list(disc_strata(7, "2", one_prong=True))), 12)
        for config in [*disc_strata(5, "2"), *closed_strata(3, "2")]:
            s = config["surface"]
            self.assertEqual(sum(p-2 for p in s["inner_prongs"]+s["outer_prongs"]), 4*s["genus"]-4)

    def test_state_budget_does_not_claim_completion(self):
        with self.assertRaises(StateBudgetExceeded):
            Automaton.construct(disc_config(limits={"state_decorations": 0}))

    def test_rooted_canonicalization_agrees_with_exhaustive_symmetry_orbits(self):
        from ttauto_pipeline.combinatorics import compositions, matchings
        from ttauto_pipeline.ribbon_states import graph_from_fans, symmetry_orders, canonical, rooted_transports
        inner = [2, 2]
        for colors in (["filled", "filled"], ["filled", "puncture"]):
            correspondence = {}
            orders = list(symmetry_orders(inner, colors))
            for degrees in compositions(6, 4):
                for pairing in matchings(tuple(range(6))):
                    slots = [0]*6
                    for e, (a, b) in enumerate(pairing, 1):
                        slots[a], slots[b] = e, -e
                    fans, at = [], 0
                    for degree in degrees:
                        fans.append(slots[at:at+degree]); at += degree
                    if graph_from_fans(inner, fans) is None:
                        continue
                    exhaustive = canonical(fans, orders)
                    rooted = min(rooted_transports(inner, colors, fans))[0]
                    correspondence.setdefault(exhaustive, set()).add(rooted)
            self.assertTrue(correspondence)
            self.assertTrue(all(len(values) == 1 for values in correspondence.values()))
            self.assertEqual(len({next(iter(values)) for values in correspondence.values()}), len(correspondence))


class PipelineTests(unittest.TestCase):
    def test_all_pruning_switches_preserve_complete_small_survivors(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline = disc_config(pruning={"norm": False})
            run(baseline, root/"off")
            wanted = survivors(records(root/"off"))
            self.assertEqual(len(wanted), 2)
            for name in ("norm", "cycle", "homology", "mod2", "lefschetz", "perron", "all"):
                switches = {key: name == "all" or key == name for key in ("norm", "cycle", "homology", "mod2", "lefschetz", "perron")}
                report = run(disc_config(pruning=switches), root/name)
                self.assertTrue(report["enumeration_completed"])
                self.assertFalse(report["minimum_proved"])
                self.assertEqual(survivors(records(root/name)), wanted, name)
            self.assertTrue(verify(root/"all", replay_search=True)["verified"])

    def test_resume_is_exact_and_changed_policy_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = disc_config(limits={"events": 9})
            report = run(config, root/"resumed", checkpoint_every=1)
            self.assertEqual(report["status"], "paused_budget")
            self.assertFalse(report["enumeration_completed"])
            config["limits"]["events"] = 10000
            resumed = run(config, root/"resumed", resume=True)
            direct = run(config, root/"direct")
            self.assertEqual(resumed["counts"], direct["counts"])
            self.assertEqual(records(root/"resumed"), records(root/"direct"))
            old = (root/"resumed"/"result.json").read_bytes()
            config["pruning"]["cycle"] = True
            with self.assertRaisesRegex(ValueError, "resume refused"):
                run(config, root/"resumed", resume=True)
            self.assertEqual((root/"resumed"/"result.json").read_bytes(), old)

    def test_geometry_failure_rolls_back_and_is_retryable(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch("ttauto_pipeline.engine.classify", side_effect=TimeoutError("control")):
                with self.assertRaises(TimeoutError):
                    run(disc_config(), folder, checkpoint_every=1)
            report = json.loads((Path(folder)/"result.json").read_text())
            self.assertEqual(report["status"], "failed")
            self.assertFalse(report["enumeration_completed"])
            self.assertNotIn("attained", report["counts"].get("geometry", {}))
            report = run(disc_config(), folder, resume=True)
            self.assertTrue(report["enumeration_completed"])
            self.assertEqual(report["counts"]["geometry"], {"attained": 2})

    def test_corrupted_record_and_automaton_are_detected(self):
        with tempfile.TemporaryDirectory() as folder:
            run(disc_config(), folder)
            with sqlite3.connect(Path(folder)/"run.sqlite3") as db:
                db.execute("UPDATE candidates SET payload='{}' WHERE id=1")
            db.close()
            with self.assertRaisesRegex(ValueError, "corrupt candidate"):
                verify(folder)
        a = Automaton.construct(disc_config())
        data = deepcopy(a.data)
        data["states"][0]["arcs"][0]["matrix"][0][0] += 1
        with self.assertRaisesRegex(ValueError, "word/matrix"):
            Automaton(data, disc_config()["surface"])

    def test_edited_minimum_claim_is_not_a_verified_report(self):
        with tempfile.TemporaryDirectory() as folder:
            report = run(disc_config(), folder)
            report["minimum_proved"] = True
            (Path(folder)/"result.json").write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, "proof claims"):
                verify(folder)

    def test_orientation_cover_and_exact_root(self):
        c = disc_config()
        a = Automaton.construct(c)
        w, v = a.compose(0, [0, 1], 0)
        event = dict(kind="candidate", start=0, path=[0, 1], terminal=0, matrix=real_matrix(w, 2))
        context = Context(a, c, event)
        self.assertEqual(context.cover["genus"], 1)
        self.assertEqual(context.cover["polynomial"], [1, -3, 1])
        self.assertEqual(sorted(context.cover["lifted_prongs"]), [2]*4)

    def test_closed_g2_geometric_attainment(self):
        config = normalize(dict(schema="ttauto-experiment-v1", surface=dict(kind="closed", genus=2, inner_prongs=[2], outer_prongs=[6]), upper="7/4"))
        a = Automaton.construct(config)
        self.assertEqual((len(a.states), sum(map(len, a.arcs))), (76, 262))
        w, v = a.compose(72, [1, 2, 1], 0)
        event = dict(kind="candidate", start=72, path=[1, 2, 1], terminal=0, matrix=real_matrix(w, 4))
        report = evaluate(a, config, event, modules(config), {})
        self.assertEqual(report["status"], "attained")
        self.assertEqual(report["polynomial"], [1, -1, -1, -1, 1])
        self.assertEqual(report["geometry"]["actual_singular_prongs"], [6])

    def test_closed_g3_actual_canonical_two_gate_junction(self):
        data = json.loads((FIXTURES/"g3_q8_witness.json").read_text())
        a = Automaton(data, data["surface"])
        for enabled in (False, True):
            config = normalize(dict(schema="ttauto-experiment-v1", surface=data["surface"], upper="21/10",
                                    pruning={name: enabled for name in ("cycle", "homology", "mod2", "lefschetz", "perron")}))
            w, v = a.compose(0, [0]*5, 0)
            event = dict(kind="candidate", start=0, path=[0]*5, terminal=0, matrix=real_matrix(w, 6))
            report = evaluate(a, config, event, modules(config), {})
            self.assertEqual(report["status"], "attained")
            self.assertEqual(report["geometry"]["actual_singular_prongs"], [10])
            self.assertEqual(len(report["geometry"]["canonical_edge_images"]), 7)
            self.assertFalse(report["geometry"]["actual_foliation_orientable"])
        config["source"] = {"kind": "automaton", "path": str(FIXTURES/"g3_q8_witness.json")}
        config["limits"] = {"depth": 5}
        with tempfile.TemporaryDirectory() as folder:
            report = run(config, folder)
            self.assertTrue(report["depth_truncated"])
            self.assertFalse(report["enumeration_completed"])
            self.assertTrue(verify(folder, replay_search=True)["verified"])


class NativeIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build = Path(os.environ.get("TTAUTO_BUILD", ROOT/"build"))
        suffix = ".exe" if os.name == "nt" else ""
        cls.export = build/"bin"/("ttauto_pipeline_export"+suffix)
        if not cls.export.is_file():
            raise unittest.SkipTest("native targets not built")

    def test_native_graph_builtin_geometry_and_independent_replay(self):
        config = disc_config(source=dict(kind="native", executable=str(self.export), stratum_index=0),
                             geometry=dict(kind="builtin"),
                             pruning={name: True for name in ("norm", "cycle", "homology", "mod2", "lefschetz", "perron")})
        with tempfile.TemporaryDirectory() as folder:
            report = run(config, folder)
            self.assertTrue(report["enumeration_completed"])
            self.assertEqual(report["counts"]["geometry"], {"attained": 2})
            self.assertEqual(report["best_attained"]["polynomial"], [1, -3, 1])
            self.assertTrue(verify(folder, replay_search=True)["verified"])

    def test_native_four_puncture_graph_matches_independent_inventory(self):
        config = normalize(dict(schema="ttauto-experiment-v1", surface=dict(kind="disc", inner_prongs=[1]*4, outer_prongs=[2]), upper="5/2"))
        independent = Automaton.construct(config)
        config["source"] = dict(kind="native", executable=str(self.export), stratum_index=0)
        native = Automaton.construct(config)
        self.assertEqual((len(native.states), sum(map(len, native.arcs))), (4, 14))
        self.assertEqual((len(independent.states), sum(map(len, independent.arcs))), (4, 14))


if __name__ == "__main__":
    unittest.main()
