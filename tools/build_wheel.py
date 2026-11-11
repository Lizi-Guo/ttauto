#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build a wheel in a clean temporary source tree, never reuse build/lib."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='ttauto-wheel-source-') as temporary:
        source = Path(temporary)
        for name in ('pyproject.toml', 'README.md', 'LICENSE', 'COPYING'):
            if (ROOT/name).is_file():
                shutil.copy2(ROOT/name, source/name)
        for package in ('ttauto_pipeline', 'ttauto_cert'):
            shutil.copytree(ROOT/'python'/package, source/'python'/package,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        wheels = source/'wheels'
        subprocess.run([sys.executable, '-m', 'pip', 'wheel', '.', '--no-build-isolation',
                        '--no-deps', '--wheel-dir', str(wheels)], cwd=source, check=True)
        for path in wheels.glob('*.whl'):
            with zipfile.ZipFile(path) as archive:
                names = archive.namelist()
                forbidden = [name for name in names if 'hall' in name.lower() or
                             name in ('ttauto_cert/bh_exact.py', 'ttauto_cert/geometry.py')]
                if forbidden:
                    raise RuntimeError(f'removed classifier modules in wheel: {forbidden}')
                expected = {f'{package}/{file.name}' for package in ('ttauto_pipeline', 'ttauto_cert')
                            for file in (ROOT/'python'/package).glob('*.py')}
                actual = {name for name in names if name.endswith('.py')}
                if actual != expected:
                    raise RuntimeError('wheel Python modules differ from current source inventory')
            shutil.copy2(path, output/path.name)
            print(f'Verified clean wheel: {output/path.name}')


if __name__ == '__main__':
    main()
