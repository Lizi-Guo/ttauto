# SPDX-License-Identifier: GPL-3.0-or-later
"""The one search pipeline: source -> prefix filters -> candidates -> geometry."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import json
import time
from ttauto_cert.exact_math import charpoly, compare_largest_real_root, compare_largest_real_roots, isolate_largest_real_root
from ttauto_cert.standard_maps import primitive
from .automaton import Automaton, FoldSource
from .config import normalize
from .geometry import classify
from .pruning import Context, modules
from .ribbon import real_matrix
from .storage import Journal, atomic_json, digest, file_digest, source_digest


def stage(counts, name, outcome):
    row = counts.setdefault("stages", {}).setdefault(name, {"tested": 0, "passed": 0, "rejected": 0, "skipped": 0})
    row["tested"] += 1
    row[{"pass": "passed", "reject": "rejected", "skip": "skipped"}[outcome]] += 1


def evaluate(automaton, config, candidate, plugins, counts):
    record = {**candidate, "decisions": []}
    context = Context(automaton, config, record)

    def decision(name, outcome, reason, evidence=None):
        stage(counts, name, outcome)
        record["decisions"].append({"module": name, "outcome": outcome, "reason": reason, "evidence": evidence})
        if outcome == "reject":
            record.update(status="rejected", rejected_by=name)
            return True
        return False

    for module in plugins:
        if module.stage == "closed_matrix":
            d = module.check(context)
            if decision(module.name, d.outcome, d.reason, d.evidence):
                return record
    if decision("primitive", "pass" if primitive(record["matrix"]) else "reject", "strong_connectivity_and_cycle_gcd"):
        return record
    polynomial = record["polynomial"] = charpoly(record["matrix"])
    within = compare_largest_real_root(polynomial, 1) > 0 and compare_largest_real_root(polynomial, Fraction(config["upper"])) <= 0
    if decision("growth", "pass" if within else "reject", "exact_1_lt_Perron_le_upper"):
        return record
    if config.get("target_polynomial") is not None:
        matched = compare_largest_real_roots(polynomial, config["target_polynomial"]) == 0
        if decision("target", "pass" if matched else "reject", "exact_largest_root_equality"):
            return record
    for module in plugins:
        if module.stage == "signed_map":
            d = module.check(context)
            if decision(module.name, d.outcome, d.reason, d.evidence):
                return record
    words, vertices = context.mapping
    if real_matrix(words, automaton.dimension) != record["matrix"]:
        raise ValueError("composed signed word disagrees with matrix product")
    record["edge_images"], record["vertex_images"] = words, vertices
    geometry = classify(context)
    record["geometry"] = geometry
    record["status"] = geometry["status"]
    counts.setdefault("geometry", {}).setdefault(geometry["status"], 0)
    counts["geometry"][geometry["status"]] += 1
    if record["status"] == "attained":
        left, right = isolate_largest_real_root(polynomial, bits=64)
        record["root_interval"] = [str(left), str(right)]
        record["display_value"] = float((left+right)/2)
    return record


def identity(config, automaton, sources):
    bound_config = deepcopy(config)
    bound_config["limits"].pop("events", None)
    bound_config["limits"].pop("state_decorations", None)
    if bound_config["source"].get("path"):
        bound_config["source"]["path"] = "content-bound-automaton"
    if bound_config["source"].get("executable"):
        bound_config["source"]["executable"] = file_digest(bound_config["source"]["executable"])
    return {"schema": "ttauto-run-identity-v1", "config": bound_config,
            "automaton_sha256": digest(automaton.data), "python_sources_sha256": sources}


def result(config, automaton, source, checkpoint, status, run_identity, error=None):
    counts = checkpoint["counts"]
    completed = status == "exhausted" and not source.depth_truncated
    report = {"schema": "ttauto-result-v1", "name": config.get("name", "experiment"),
              "status": status, "enumeration_completed": completed,
              "depth_truncated": source.depth_truncated, "norm_bound": source.bound,
              "depth_bound": source.depth, "theoretical_depth_bound": source.theoretical_depth,
              "surface": config["surface"], "upper": config["upper"],
              "states": len(automaton.states), "arcs": sum(map(len, automaton.arcs)),
              "terminal_symmetries": sum(map(len, automaton.terminals)),
              "counts": counts, "best_attained": checkpoint.get("best"),
              "minimum_proved": False,
              "proof_scope": "bounded paths in the supplied automaton, with every terminal symmetry; whole-surface normal-form coverage is a separate theorem obligation",
              "unresolved_candidates": counts.get("geometry", {}).get("unresolved", 0),
              "run_identity": run_identity, "checkpoint_sha256": digest(checkpoint)}
    if error:
        report["error"] = error
    return report


def run(config, output, *, resume=False, progress=None, checkpoint_every=100):
    config = normalize(config)
    folder = Path(output).resolve()
    if not resume and folder.exists() and any(folder.iterdir()):
        raise ValueError("output directory is not empty; select a new directory or use --resume")
    folder.mkdir(parents=True, exist_ok=True)
    before = source_digest()
    backend = None
    graph_path = folder/"automaton.json"
    if config['schema'] == 'ttauto-proof-v1':
        from .proof_backend import ProofBackend
        backend = ProofBackend(config, folder, resume, before)
        automaton = None
    elif resume:
        if not graph_path.is_file():
            raise ValueError("resume requires automaton.json")
        automaton = Automaton(json.loads(graph_path.read_text(encoding="utf-8")), config["surface"])
        if config["source"]["kind"] == "automaton":
            supplied = json.loads(Path(config["source"]["path"]).read_text(encoding="utf-8"))
            if digest(supplied) != digest(automaton.data):
                raise ValueError("resume refused: supplied automaton changed")
    else:
        automaton = Automaton.construct(config)
        atomic_json(graph_path, automaton.data)
    run_identity = backend.identity if backend else identity(config, automaton, before)
    journal = Journal(folder/"run.sqlite3", run_identity, resume=resume)
    try:
        checkpoint = journal.get("checkpoint") or {"cursor": None, "counts": {"events": 0, "prefixes": 0, "candidates": 0}, "best": None}
        source = backend.source(checkpoint['cursor']) if backend else FoldSource(automaton, config, checkpoint["cursor"])
        plugins = [] if backend else modules(config)
        status = "running"
        limit = config["limits"].get("events")
        tick = time.monotonic()
        # Ensure a valid restart point exists even before the first candidate.
        checkpoint["cursor"] = source.cursor
        checkpoint["status"] = "running"
        checkpoint.pop("error", None)
        journal.checkpoint(checkpoint)
        atomic_json(folder/"config.json", config)

        def save(state):
            checkpoint["cursor"] = source.cursor
            checkpoint["status"] = state
            journal.checkpoint(checkpoint)
            report = backend.report(source, checkpoint, state) if backend else result(config, automaton, source, checkpoint, state, run_identity)
            atomic_json(folder/"result.json", report)
            return report

        try:
            while limit is None or checkpoint["counts"]["events"] < limit:
                event = source.next_event()
                if event is None:
                    status = "exhausted"
                    break
                counts = checkpoint["counts"]
                counts["events"] += 1
                if event["kind"] == "prefix":
                    counts["prefixes"] += 1
                    if backend:
                        backend.record_prefix(event, counts)
                    if not backend and config["pruning"]["norm"]:
                        rejected = sum(map(sum, event["matrix"])) > source.bound
                        stage(counts, "norm", "reject" if rejected else "pass")
                        if rejected:
                            source.reject_prefix()
                else:
                    counts["candidates"] += 1
                    record = backend.evaluate(event, counts) if backend else evaluate(automaton, config, event, plugins, counts)
                    journal.add_candidate(counts["candidates"], record)
                    best = checkpoint["best"]
                    if record["status"] == "attained" and (best is None or compare_largest_real_roots(record["polynomial"], best["polynomial"]) < 0):
                        checkpoint["best"] = {key: record[key] for key in ("start", "path", "terminal", "polynomial", "root_interval", "display_value")}
                if counts["events"] % checkpoint_every == 0:
                    report = save("running")
                    if progress and time.monotonic()-tick >= 5:
                        progress(report)
                        tick = time.monotonic()
            else:
                status = "paused_budget"
            if backend:
                backend.validate_runtime()
            if source_digest() != before or (not backend and identity(config, automaton, before) != run_identity):
                raise ValueError("code or backend changed during the run")
            return save(status)
        except BaseException as error:
            # Roll back the current batch and retry it on resume. A timeout,
            # interrupt or arithmetic failure is never a mathematical exclusion.
            journal.connection.rollback()
            checkpoint = journal.get("checkpoint")
            source = backend.source(checkpoint['cursor']) if backend else FoldSource(automaton, config, checkpoint["cursor"])
            checkpoint["status"] = "failed"
            checkpoint["error"] = f"{type(error).__name__}: {error}"
            journal.checkpoint(checkpoint)
            report = (backend.report(source, checkpoint, 'failed', checkpoint['error']) if backend else
                      result(config, automaton, source, checkpoint, "failed", run_identity, checkpoint["error"]))
            atomic_json(folder/"result.json", report)
            raise
    finally:
        journal.close()
