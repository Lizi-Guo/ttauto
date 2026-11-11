# SPDX-License-Identifier: GPL-3.0-or-later
"""Independent attainment, scope boundaries and geometric evidence replay."""
from copy import deepcopy
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from ttauto_cert.exact_math import charpoly
from ttauto_pipeline.automaton import Automaton
from ttauto_pipeline.config import normalize
from ttauto_pipeline.engine import run
from ttauto_pipeline.geometry import classify, pillowcase
from ttauto_pipeline.pruning import Context
from ttauto_pipeline.ribbon import real_matrix
from ttauto_pipeline.storage import canonical, digest
from ttauto_pipeline.verify import verify


def context():
    config = normalize(dict(schema="ttauto-experiment-v1", surface=dict(
        kind="disc", inner_prongs=[1]*3, outer_prongs=[1]), upper="8/3"))
    automaton = Automaton.construct(config)
    words, _ = automaton.compose(0, [0, 1], 0)
    matrix = real_matrix(words, 2)
    return Context(automaton, config, dict(start=0, path=[0, 1], terminal=0,
                                         matrix=matrix, polynomial=charpoly(matrix)))


class BuiltinGeometryTests(unittest.TestCase):
    def test_complete_disc_run_without_any_external_process(self):
        with tempfile.TemporaryDirectory() as folder, patch(
                "subprocess.Popen", side_effect=AssertionError("external process forbidden")):
            report = run(context().config, folder)
            self.assertEqual(report["counts"]["geometry"], {"attained": 2})
            self.assertEqual(report["best_attained"]["polynomial"], [1, -3, 1])
            self.assertTrue(verify(folder, replay_search=True)["verified"])

    def test_certificate_is_for_the_actual_signed_cover_action(self):
        c = context()
        result = classify(c)
        self.assertEqual(result["cover"], c.cover)
        self.assertEqual(result["trace"], 3)
        self.assertEqual(result["determinant"], 1)
        self.assertEqual(result["actual_puncture_prongs"], [1]*3)
        self.assertFalse(result["actual_foliation_orientable"])

    def test_audit_does_not_claim_attainment(self):
        c = context()
        c.config["geometry"] = {"kind": "audit"}
        with tempfile.TemporaryDirectory() as folder:
            report = run(c.config, folder)
            self.assertEqual(report["counts"]["geometry"], {"unresolved": 2})
            self.assertIsNone(report["best_attained"])
            self.assertTrue(verify(folder, replay_search=True)["verified"])

    def test_fourth_disc_puncture_is_outside_the_certificate(self):
        c = context()
        # Even a hyperbolic torus matrix must not certify an extra puncture.
        c._cover = deepcopy(c.cover)
        c.config = normalize(dict(schema="ttauto-experiment-v1", surface=dict(
            kind="disc", inner_prongs=[1]*4, outer_prongs=[2]), upper="8/3"))
        self.assertEqual(classify(c)["status"], "unresolved")

    def test_elliptic_and_parabolic_cover_actions_are_not_attainers(self):
        for matrix in ([[0, -1], [1, 0]], [[1, 1], [0, 1]]):
            c = context()
            c._cover = {**c.cover, "positive_lift_matrix": matrix}
            self.assertEqual(pillowcase(c)["reason"], "nonhyperbolic_torus_action")

    def test_other_lift_has_the_same_dilatation(self):
        c = context()
        c._cover = {**c.cover, "positive_lift_matrix": [
            [-x for x in row] for row in c.cover["positive_lift_matrix"]]}
        result = pillowcase(c)
        self.assertEqual(result["status"], "attained")
        self.assertEqual(result["actual_dilatation_polynomial"], [1, -3, 1])

    def test_wrong_dilatation_and_orientability_do_not_match(self):
        c = context()
        c.candidate["polynomial"] = [1, -4, 1]
        self.assertEqual(pillowcase(c)["reason"], "actual_dilatation_mismatch")
        c = context()
        c.config["surface"]["foliation"] = "orientable"
        self.assertEqual(pillowcase(c)["reason"], "actual_stratum_mismatch")

    def test_invalid_cover_is_an_error_not_an_exclusion(self):
        for changes in ({"genus": 2}, {"lifted_prongs": [4]*4},
                        {"positive_lift_matrix": [[2, 0], [0, 2]]}):
            c = context()
            c._cover = {**c.cover, **changes}
            with self.assertRaises(ValueError):
                pillowcase(c)

    def test_external_geometry_configuration_is_rejected(self):
        for geometry in ({"kind": "hall"}, {"kind": "builtin", "executable": "unused"},
                         {"kind": "builtin", "timeout": 1}):
            config = context().config
            config["geometry"] = geometry
            with self.assertRaises(ValueError):
                normalize(config)

    def test_all_geometric_evidence_is_replayed_including_extra_reports(self):
        with tempfile.TemporaryDirectory() as folder:
            run(context().config, folder)
            with sqlite3.connect(Path(folder)/"run.sqlite3") as db:
                for number, payload in db.execute("SELECT id,payload FROM candidates").fetchall():
                    record = json.loads(payload)
                    if record["status"] == "attained":
                        record["geometry"]["report"] = {"unverified": True}
                        db.execute("UPDATE candidates SET payload=?,digest=? WHERE id=?",
                                   (canonical(record), digest(record), number))
                        break
            db.close()
            with self.assertRaisesRegex(ValueError, "geometry"):
                verify(folder)


if __name__ == "__main__":
    unittest.main()
