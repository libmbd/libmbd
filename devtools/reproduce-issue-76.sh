#!/bin/bash
# Reproduce https://github.com/libmbd/libmbd/issues/76
#
# The ScaLAPACK/MPI gradient tests (e.g. mbd_rsscs_deriv_expl) crash on the
# big-endian s390x architecture with an MPICH error like:
#
#   MPIDI_POSIX_mpi_bcast_release_gather(132)..:
#   MPIDI_POSIX_mpi_release_gather_release(218): message sizes do not match
#       across processes in the collective routine: Received 0 but expected 216
#
# Root cause is NOT in libMBD. It is a big-endian bug in MPICH's ch4 POSIX
# shared-memory broadcast ("release_gather"), which mis-handles the strided
# vector datatype that ScaLAPACK (PDGETRF/PDGETRI, used by libMBD's distributed
# matrix inverse) broadcasts during the rSSCS derivative tests. See the notes in
# devtools/issue-76-notes.md for the full analysis and the run-time workaround.
#
# This script reproduces the crash inside an emulated big-endian s390x container.
# It needs Docker with QEMU/binfmt set up for foreign architectures:
#
#   docker run --privileged --rm tonistiigi/binfmt --install s390x
#
# Then, from the repository root:
#
#   devtools/reproduce-issue-76.sh            # reproduce the crash (expected fail)
#   WORKAROUND=1 devtools/reproduce-issue-76.sh   # apply the MPICH workaround (passes)
#
# A Fedora image is used because its stock MPICH is built with the ch4 device
# (ch4:ofi), which is what triggers the bug. Debian/Ubuntu MPICH uses the older
# ch3 device, which does NOT exhibit the issue, so those images cannot reproduce
# it even on s390x.
set -euo pipefail

IMAGE=${IMAGE:-fedora:40}
NRANKS=${NRANKS:-3}
TEST=${TEST:-mbd_rsscs_deriv_expl}
WORKAROUND=${WORKAROUND:-0}

repo_root=$(cd "$(dirname "$0")/.." && pwd)

cvar=""
if [ "$WORKAROUND" = 1 ]; then
    # Force the generic broadcast algorithm so MPICH bypasses the buggy
    # big-endian ch4 POSIX release_gather shared-memory path.
    cvar="MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM=mpir"
fi

exec docker run --rm --platform linux/s390x \
    -v "$repo_root":/src:ro \
    -e CVAR="$cvar" -e NRANKS="$NRANKS" -e TEST="$TEST" \
    "$IMAGE" bash -euxo pipefail -c '
        dnf install -y -q mpich-devel scalapack-mpich-devel gcc-gfortran \
            lapack-devel cmake make git diffutils
        export PATH=/usr/lib64/mpich/bin:$PATH
        export LD_LIBRARY_PATH=/usr/lib64/mpich/lib:${LD_LIBRARY_PATH:-}
        mpichversion | grep -iE "MPICH (Version|Device)"

        cp -r /src /build && cd /build && rm -rf build
        cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_SCALAPACK_MPI=ON \
            -DCMAKE_Fortran_COMPILER=mpif90 -DMPI_Fortran_COMPILER=mpif90 \
            -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
            -DSCALAPACK_LIBRARY=/usr/lib64/mpich/lib/libscalapack.so
        cmake --build build -j "$(nproc)"

        env $CVAR mpiexec -n "$NRANKS" ./build/tests/mbd_grad_tests "$TEST"
    '
