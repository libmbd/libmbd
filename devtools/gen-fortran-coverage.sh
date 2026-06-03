#!/usr/bin/env bash
# Generate .gcov coverage reports for the libMBD Fortran sources.
#
# Codecov's built-in gcov pass (and gcovr) mis-parse the cpp-preprocessed
# *.F90 sources: they report only a fraction of the executable lines and zero
# out fully exercised files (e.g. mbd_coulomb, mbd_c_api). Invoking `gcov`
# directly maps the preprocessed line table back onto the original source
# correctly, so we produce the reports here and hand them to Codecov with its
# own gcov plugin disabled (plugins: noop).
#
# Usage: gen-fortran-coverage.sh [OUTDIR]
#   BLDDIR  build directory containing the .gcda data (default: build)
#   GCOV    gcov executable to use (default: derived from $FC / gfortran)
set -euo pipefail

BLDDIR=${BLDDIR:-build}
outdir=${1:-coverage-fortran}

# A gcov that does not match the compiler silently emits empty data
# ("version 'B13*', prefer 'A13*'"), so derive it from the Fortran compiler:
# gfortran -> gcov, gfortran-13 -> gcov-13, <triple>-gfortran -> <triple>-gcov.
if [[ -z "${GCOV:-}" ]]; then
    fc=$(command -v "${FC:-gfortran}" 2>/dev/null || true)
    if [[ -n "$fc" ]]; then
        GCOV=$(printf '%s' "$fc" | sed -E 's#gfortran(-[0-9]+)?$#gcov\1#')
    else
        GCOV=gcov
    fi
fi
command -v "$GCOV" >/dev/null 2>&1 || { echo "gcov executable '$GCOV' not found" >&2; exit 1; }

# Only the library objects matter; the test-executable .gcda cover the test
# sources, which are ignored in codecov.yml.
absbld=$(cd "$BLDDIR" && pwd)
mapfile -t gcda < <(find "$absbld" -path '*mbd.dir*' -name '*.gcda')
if [[ ${#gcda[@]} -eq 0 ]]; then
    echo "no .gcda coverage data found under $BLDDIR" >&2
    exit 1
fi

# gcov writes .gcov files into the current directory; -p keeps the full path in
# the name and -b adds branch data. Collect everything in $outdir.
mkdir -p "$outdir"
( cd "$outdir" && "$GCOV" -pb "${gcda[@]}" ) >/dev/null
echo "generated $(find "$outdir" -name '*.gcov' | wc -l) .gcov file(s) in $outdir using $GCOV"
