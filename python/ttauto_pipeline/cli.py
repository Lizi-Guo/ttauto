# SPDX-License-Identifier: GPL-3.0-or-later
"""One CLI for every surface, pruning combination and replay."""
import argparse
import json
from pathlib import Path
import sys
from . import __version__
from .automaton import Automaton
from .config import load
from .engine import run
from .pruning import REGISTRY
from .storage import atomic_json
from .verify import verify
from .catalog import disc_strata, closed_strata
from .campaign import campaign
from .proof_config import MODULES as PROOF_MODULES
from .proof_pruning import CONTRACTS


def main(argv=None):
    parser = argparse.ArgumentParser(prog="ttauto", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "build", "validate"):
        p = commands.add_parser(name)
        p.add_argument("config", type=Path)
        if name != "validate":
            p.add_argument("--output", type=Path, required=True)
        if name == "run":
            p.add_argument("--resume", action="store_true")
            p.add_argument("--events", type=int, help="total event budget, including already committed events")
            choices = sorted({'norm', *REGISTRY, *PROOF_MODULES})
            p.add_argument("--enable", action="append", choices=choices, default=[])
            p.add_argument("--disable", action="append", choices=choices, default=[])
    p = commands.add_parser("verify")
    p.add_argument("run", type=Path)
    p.add_argument("--replay-search", action="store_true")
    p = commands.add_parser("summary")
    p.add_argument("run", type=Path)
    commands.add_parser("modules")
    p = commands.add_parser('import-catalog')
    p.add_argument('input', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--field', help='explicit dotted JSON field containing coefficient lists, e.g. survivors')
    p = commands.add_parser("catalog")
    p.add_argument("surface", choices=("disc", "closed"))
    p.add_argument("size", type=int, help="puncture count or genus")
    p.add_argument("--upper", required=True)
    p.add_argument("--one-prong", action="store_true")
    p.add_argument("--auxiliary-points", type=int, default=1)
    p.add_argument("--output", type=Path, required=True)
    p = commands.add_parser("campaign")
    p.add_argument("configs", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--events", type=int)
    args = parser.parse_args(argv)
    try:
        if args.command == "modules":
            result = [{"name": "norm", "stage": "prefix", "contract": "permutation-plus-one folds; exact total-entry bound"}]
            result += [{"name": name, "stage": cls.stage, "contract": cls.contract} for name, cls in REGISTRY.items()]
            result += [{'name': name, 'stage': 'spectral_proof', 'contract': contract} for name, contract in CONTRACTS.items()]
        elif args.command == 'import-catalog':
            from .import_catalog import import_catalog
            result = import_catalog(args.input, args.output, args.field)
        elif args.command == "campaign":
            result = campaign(args.configs, args.output, args.resume, args.events)
        elif args.command == "catalog":
            if args.output.exists() and any(args.output.iterdir()):
                raise ValueError("catalog output must be empty")
            args.output.mkdir(parents=True, exist_ok=True)
            specs = disc_strata(args.size, args.upper, args.one_prong) if args.surface == "disc" else closed_strata(args.size, args.upper, args.auxiliary_points)
            count = 0
            for count, spec in enumerate(specs, 1):
                atomic_json(args.output/f"{count:05d}.json", spec)
            result = {"configurations": count, "scope": "formal prong types; existence and whole-surface coverage are not asserted"}
        elif args.command == "verify":
            result = verify(args.run, args.replay_search)
        elif args.command == "summary":
            result = json.loads((args.run/"result.json").read_text(encoding="utf-8"))
        else:
            config = load(args.config)
            if args.command == "validate":
                result = config
            elif args.command == "build":
                if args.output.exists():
                    raise ValueError("build output already exists")
                args.output.parent.mkdir(parents=True, exist_ok=True)
                if config['schema'] == 'ttauto-proof-v1':
                    import tempfile
                    from .proof_backend import ProofBackend
                    from .storage import source_digest
                    with tempfile.TemporaryDirectory() as temporary:
                        backend = ProofBackend(config, temporary, False, source_digest())
                        atomic_json(args.output, backend.space)
                        result = {'strata': len(backend.space['cases']), 'automaton_required': False, 'output': str(args.output)}
                else:
                    automaton = Automaton.construct(config)
                    atomic_json(args.output, automaton.data)
                    result = {"states": len(automaton.states), "arcs": sum(map(len, automaton.arcs)), "output": str(args.output)}
            else:
                if set(args.enable) & set(args.disable):
                    raise ValueError("a module cannot be both enabled and disabled")
                for name in args.enable:
                    config["pruning"][name] = True
                for name in args.disable:
                    config["pruning"][name] = False
                if args.events is not None:
                    config["limits"]["events"] = args.events
                def progress(report):
                    print(json.dumps({"status": report["status"], "events": report["counts"]["events"],
                                      "candidates": report["counts"]["candidates"]}), file=sys.stderr, flush=True)
                result = run(config, args.output, resume=args.resume, progress=progress)
        print(json.dumps(result, indent=2))
        return 0
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        print(f"ttauto: {type(error).__name__}: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
