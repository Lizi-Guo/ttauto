# Preparing a GitHub release

This directory contains the algorithms and numerical libraries as ordinary
source files. It has no submodules and needs no upstream checkout to build.
Publish it to a repository you own. Original project URLs record attribution;
they are not dependency-fetch endpoints. No upload is performed by the tools.

Before a release:

1. Build and run CTest, including the Python pipeline and native integration tests.
2. Run `experiments/proofs/disc3-minimum.json` and
   `experiments/proofs/genus3-bulk21.json`, then replay their complete searches.
   Run a budgeted native large-disc example and the explicit-map regressions.
3. Inspect `git diff` and `git status`. Keep historical campaign folders, full
   search databases, logs, user paths, credentials, and binaries out of source control.
4. Run `python tools/release.py --output ../ttauto-source.zip`. The allowlisted
   archive includes sources, documentation, small fixtures, dependency licenses
   and a checksum manifest. A ZIP consumer can build with CMake and a compiler
   without Git, a submodule checkout or network access. For Python execution
   without installation, set `PYTHONPATH=python` and use `python -m ttauto_pipeline`.
   Build Python wheels with `python tools/build_wheel.py --output dist`. It uses
   a fresh source tree and compares the wheel's modules with the current source
   inventory, preventing deleted classifiers in an old `build/lib` cache from
   reappearing in a release.
5. Create a repository under the intended GitHub account, add that repository as
   a remote, commit the reviewed source changes and push explicitly.

The source archive excludes `runs/` and raw historical campaign archives; the
documentation index `docs/HISTORICAL_RESULTS.json` is included. Release descriptions
must distinguish the new implementation's verified examples from earlier
distributed computations. Do not claim the full D10/D12/D14/D18/D22/D26 or all
genus-three minima were rerun simply because configurations now exist.

The original C++ search remains for compatibility and comparison. The supported
experiment entrypoint is `ttauto`; use it for new shared-pipeline reports.

Source identities changed in 0.1.1, 0.1.2 and 0.2.0. Keep earlier run directories as
historical records and use new output directories; do not resume them under
changed code. Version 0.1.2 has no external geometry backend. Old Hall-based
reports remain historical evidence and are not accepted as new certificates.
Version 0.2.0 adds spectral and verified-residual sources to the shared engine,
ports the old large-disc/genus-three arithmetic kernels, and fixes actual
Lefschetz lift actions to use JA. Spectral lower-bound certificates and exact
minimum certificates have different requirements; preserve those distinctions
in release notes. The source examples' target polynomials are input questions,
not tabulated answers substituted for computation.
