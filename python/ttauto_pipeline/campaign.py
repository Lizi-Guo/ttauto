# SPDX-License-Identifier: GPL-3.0-or-later
"""Execute any collection of strata through the same run engine."""
from pathlib import Path
from .config import load
from .engine import run
from .storage import atomic_json


def campaign(config_folder, output, resume=False, events=None, progress=None):
    paths = sorted(Path(config_folder).glob("*.json"))
    if not paths:
        raise ValueError("no JSON experiment configurations found")
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()) and not resume:
        raise ValueError("campaign output exists; use --resume or a new directory")
    output.mkdir(parents=True, exist_ok=True)
    report = {"schema": "ttauto-campaign-v1", "experiments": [], "completed": False,
              "minimum_proved": False}
    for path in paths:
        config = load(path)
        if events is not None:
            config["limits"]["events"] = events
        destination = output/path.stem
        try:
            result = run(config, destination, resume=resume and (destination/"run.sqlite3").exists(), progress=progress)
        except Exception as error:
            report["error"] = {"configuration": path.name, "detail": f"{type(error).__name__}: {error}"}
            atomic_json(output/"campaign.json", report)
            raise
        report["experiments"].append({"configuration": path.name, "result": f"{path.stem}/result.json",
                                      "status": result["status"], "enumeration_completed": result["enumeration_completed"],
                                      "lower_bound_proved": result.get('lower_bound_proved', False),
                                      "minimum_proved": result['minimum_proved'],
                                      "counts": result["counts"], "best_attained": result["best_attained"]})
        atomic_json(output/"campaign.json", report)
    report["completed"] = all(row["enumeration_completed"] for row in report["experiments"])
    report['all_lower_bounds_proved'] = all(row['lower_bound_proved'] for row in report['experiments'])
    report['all_minima_proved'] = all(row['minimum_proved'] for row in report['experiments'])
    atomic_json(output/"campaign.json", report)
    return report
