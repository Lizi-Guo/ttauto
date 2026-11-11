# SPDX-License-Identifier: GPL-3.0-or-later
"""Replay run records; optionally re-enumerate the entire bounded search."""
import json
from pathlib import Path
import sqlite3
import tempfile
from .automaton import Automaton, FoldSource
from .config import normalize
from .engine import evaluate, identity, run, result as make_result
from .pruning import modules
from .storage import Journal, digest, source_digest


def verify(folder, replay_search=False):
    folder = Path(folder).resolve()
    config = normalize(json.loads((folder/"config.json").read_text(encoding="utf-8")))
    if config['schema'] == 'ttauto-proof-v1':
        from .proof_verify import verify as verify_proof
        return verify_proof(folder, replay_search)
    saved = json.loads((folder/"result.json").read_text(encoding="utf-8"))
    automaton = Automaton(json.loads((folder/"automaton.json").read_text(encoding="utf-8")), config["surface"])
    run_identity = identity(config, automaton, source_digest())
    if saved["run_identity"] != run_identity:
        raise ValueError("verification refused: source, configuration, backend or automaton changed")
    journal = Journal(folder/"run.sqlite3", run_identity, resume=True)
    try:
        checkpoint = journal.get("checkpoint")
        if saved["checkpoint_sha256"] != digest(checkpoint) or saved["counts"] != checkpoint["counts"]:
            raise ValueError("report/checkpoint mismatch")
        source = FoldSource(automaton, config, checkpoint["cursor"])
        expected = make_result(config, automaton, source, checkpoint, checkpoint["status"], run_identity, checkpoint.get("error"))
        if saved != expected:
            raise ValueError("report status, proof claims or metadata differ from the committed checkpoint")
        counts, total = {}, 0
        plugins = modules(config)
        for number, record in journal.records():
            total += 1
            if number != total:
                raise ValueError("candidate sequence is not contiguous")
            words, _ = automaton.compose(record["start"], record["path"], record["terminal"])
            from .ribbon import real_matrix
            if real_matrix(words, automaton.dimension) != record["matrix"]:
                raise ValueError("candidate matrix differs from signed path replay")
            candidate = {k: record[k] for k in ("kind", "start", "path", "terminal", "matrix")}
            replay = evaluate(automaton, config, candidate, plugins, counts)
            for key in ("status", "rejected_by", "polynomial", "edge_images", "vertex_images", "decisions", "root_interval", "geometry", "display_value"):
                if digest(replay.get(key)) != digest(record.get(key)):
                    raise ValueError(f"candidate {number}: replay differs in {key}")
        if total != saved["counts"]["candidates"]:
            raise ValueError("candidate count differs from journal")
        for name, row in counts.get("stages", {}).items():
            if row != saved["counts"]["stages"][name]:
                raise ValueError(f"stage count mismatch: {name}")
        if counts.get("geometry", {}) != saved["counts"].get("geometry", {}):
            raise ValueError("geometry count mismatch")
    finally:
        journal.close()
    if replay_search:
        with tempfile.TemporaryDirectory(prefix="ttauto-replay-") as directory:
            replay = run(config, directory)
        if replay["counts"] != saved["counts"] or replay["best_attained"] != saved["best_attained"] or replay["status"] != saved["status"]:
            raise ValueError("full search replay differs")
    return {"verified": True, "candidates_replayed": total, "full_search_replayed": replay_search,
            "minimum_proved": False}
