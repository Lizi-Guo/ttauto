#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build a source-only, allowlisted, reproducible release ZIP. Never uploads."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = ("include", "lib/traintracks", "examples", "tests", "testsuite", "python",
               "cmake", "native", "docs", "doc", "experiments", "tools", ".github",
               "extern/jlt")
ROOT_FILES = ("CMakeLists.txt", "pyproject.toml", "README.md", "README", "LICENSE", "COPYING",
              ".gitignore", "dependencies.lock.json", "LOCAL_PROVENANCE.json")
SUFFIXES = {".py", ".cpp", ".c", ".cc", ".h", ".hpp", ".tcc", ".cmake", ".json", ".md",
            ".txt", ".tex", ".bib", ".sh", ".yml", ".yaml", ".toml", ".in", ".inc"}
NAMES = {"LICENSE", "COPYING", "README", "CMakeLists.txt", "Makefile", ".gitignore", "SConscript", "SConstruct"}
EXCLUDE = {".git", "__pycache__", "runs", "node_modules", ".venv", "dist"}


def files():
    selected = {ROOT/name for name in ROOT_FILES if (ROOT/name).is_file()}
    for directory in DIRECTORIES:
        for path in (ROOT/directory).rglob("*"):
            relative = path.relative_to(ROOT)
            if any(part in EXCLUDE or part == "build" or part.startswith("build-") or part.endswith(".egg-info")
                   for part in relative.parts[:-1]):
                continue
            if path.is_file() and (path.suffix in SUFFIXES or path.name in NAMES):
                if path.is_symlink() or not path.resolve().is_relative_to(ROOT):
                    raise ValueError(f"unsafe release path: {relative}")
                selected.add(path)
    return sorted(selected, key=lambda path: path.relative_to(ROOT).as_posix())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    contents = {path.relative_to(ROOT).as_posix(): path.read_bytes() for path in files()}
    manifest = {name: hashlib.sha256(data).hexdigest() for name, data in contents.items()}
    contents["RELEASE_MANIFEST.json"] = (json.dumps(manifest, indent=2)+"\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(contents.items()):
            info = zipfile.ZipInfo("ttauto/"+name, (2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100755 if name.endswith(".sh") else 0o100644) << 16
            archive.writestr(info, data)
    checksum = hashlib.sha256(args.output.read_bytes()).hexdigest()
    print(json.dumps({"archive": str(args.output.resolve()), "files": len(contents), "sha256": checksum}))


if __name__ == "__main__":
    main()
