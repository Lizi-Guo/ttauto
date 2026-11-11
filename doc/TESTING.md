# How to test the ttauto code

## Build first (CMake)

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

This compiles:

- all `examples/*.cpp` into `build/examples/`
- all `tests/*.cpp` into `build/tests/`
- all `testsuite/**/*.cpp` into `build/testsuite/`
- the library into `build/lib/`, and the two pipeline helpers
  (`ttauto_pipeline_export`, `ttauto_trace_prefix`) into `build/bin/`

Nothing is written into the source tree.  Paths such as `examples/ttauto`
below name the program; the binary is `build/examples/ttauto`.

`ctest` runs deterministic fast tests by default: the testsuite programs
and the Python pipeline tests (`pipeline_python`, `python/tests/test_*.py`).

Testsuite programs are compiled with `-UNDEBUG` (see `ttauto_add_testsuite`
in `CMakeLists.txt`), so their `assert()` calls stay live in the default
Release build.  Before this was added every assert-based test passed
vacuously.  New tests should prefer `CHECK(cond)` / `CHECK_MSG(cond, msg)`
from `testsuite/check.hpp`, which print the failing expression and its
location and exit nonzero.

```bash
ctest --test-dir build --output-on-failure
```

To include slow integration tests, configure with:

```bash
cmake -S . -B build -DTTAUTO_ENABLE_SLOW_TESTS=ON
ctest --test-dir build --output-on-failure
```

To run only slow tests:

```bash
ctest --test-dir build --output-on-failure -L slow
```

Fail-fast paths are normally left untested, since they guard invariants a
caller cannot reach.  The coding reader is the exception: its errors come
from user input, so `testsuite/traintracks/test_coding_io.cpp` runs each
bad coding in a forked child and checks the exit status and the message.
That is the only test of its kind; it is skipped where `fork()` is
unavailable.

The `tests/` folder remains available for direct-run feature and exploratory
programs.

For header/component coverage mapping, see `testsuite/COVERAGE.md`.

## Programs that can be used as tests but take a long time to run

`examples/ttauto_count`   (2 minutes)

`examples/ttauto_labels` belonged here at about 100 minutes until it began
capping its path length.  The cost was the bound `check_norms()` derives
from the dilatation window, 734 folds on that stratum; searching to it is
what makes the search complete for the window, so it was the price of
completeness rather than waste.  The classes found are unchanged at every
cap from 10 to 24, so the example caps at 12 and finishes in a fiftieth of
a second.  That cap is measured, not proved -- see
`devel/iss016/labelled_automaton.md`.

## Braids from folding paths

`examples/ttbraid <punctures> <stratum> <vertex> <branch>...` reads the
braid off one closed folding path and says whether it checks out; the
check is that the braid's growth under the Dynnikov action equals the
Perron root of the path's transition matrix.  `examples/ttbraid_strata`
does the same for the minimiser of every stratum up to seven punctures
and compares each against the published table, writing the markdown that
`examples/ttauto_strata_braids.md` holds; it takes about two and a half
seconds and is regenerated in the testsuite.

## Measuring line coverage

Use a separate Debug build with gcov instrumentation; never the Release
`build/` directory.  Every output lands under its own build directory, so
the coverage build does not overwrite the Release library or binaries and
`build/` needs no repair afterwards.

```bash
cmake -S . -B build-cov -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_FLAGS=--coverage -DCMAKE_CXX_FLAGS=--coverage \
      -DCMAKE_EXE_LINKER_FLAGS=--coverage
cmake --build build-cov -j --target $(cmake --build build-cov --target help \
      | grep -o 'testsuite_[a-z_0-9]*' | sort -u)
ctest --test-dir build-cov -LE slow
gcov --json-format -o build-cov/CMakeFiles/ttauto.dir/lib/traintracks \
      build-cov/CMakeFiles/ttauto.dir/lib/traintracks/*.gcda
```

Building only the `testsuite_*` targets is enough for the measurement and
keeps the coverage build short.  Remove it afterwards with

```bash
rm -rf build-cov
```

Header-only code (`include/ttauto/**`)
is covered through the testsuite objects
(`build-cov/CMakeFiles/testsuite_*.dir/`); merge the per-object JSON by
line to get a per-file figure.  The per-file table as of 2026-09-19 is at
the top of `testsuite/COVERAGE.md`.

## Large automata

Building the automaton for a large stratum used to segfault, because
construction recursed once per vertex (issue #14).  It now uses a
worklist and indexes vertices by their coding, so the nine-puncture
stratum `1. 1. 1. 1. 1. 1. 1. 1. 1. 3 3 (5)` builds on the default
stack: 71253 vertices in about 165 s and 1.4 GB.  Memory, not the
stack, is now what limits how far this can be pushed.

`testsuite/ttauto/test_ttfoldgraph_build.cpp` guards the property by
building a seven-puncture stratum on a thread with a 256 KB stack; it
is skipped where pthreads are unavailable.

## Gate-test census

`tests/ttauto_gate_census [nmin [nmax]]` (default 3 6; n=7 takes about
9 s) repeats the strata scan with the Bestvina-Handel gate test on and lists
every class rejected by it, per stratum, with dilatation, characteristic
polynomial and representative paths.  It also counts, among the stored
representatives of the accepted classes, those whose gate partition at an
unpunctured multigon is finer than its prongs (zero on every stratum up to
n=7 as of 2026-09-19).  Use it after any change to the fold map or the gate
test; the expected rejections are recorded in
`devel/iss002/ISSUE2_IMPLEMENTATION_PLAN.md` (Step 4 outcome).

## Programs that should be used for testing

`tests/test_permplus1`
`tests/test_traintrack`
`tests/test_folding_path`
`tests/test_badwords`
`examples/ttbraid`
`examples/ttauto_min_example`
`examples/ttauto_torus`

The program `examples/ttauto` has interactive input but the defaults can be
used:
```
build/examples/ttauto << EOF
EOF
```
This produces a file `pA_n=5_1_1_inv.m` in the current directory; delete it
before committing.
