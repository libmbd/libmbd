# Issue #76: ScaLAPACK/MPI gradient tests crash on big-endian s390x

Reference: <https://github.com/libmbd/libmbd/issues/76>

## Symptom

Built with `-DENABLE_SCALAPACK_MPI=ON` and run under MPICH on s390x
(big-endian), the rSSCS derivative gradient tests abort, e.g.
`mbd_rsscs_deriv_expl` with `BLACS grid: 1 x 3`:

```
MPIDI_POSIX_mpi_bcast_release_gather(132)..:
MPIDI_POSIX_mpi_release_gather_release(218): message sizes do not match
    across processes in the collective routine: Received 0 but expected 216
```

The same build and tests pass on little-endian (x86-64, aarch64).

## Root cause (not in libMBD)

The failure is a **big-endian bug in MPICH's `ch4` device**, specifically the
POSIX shared-memory broadcast fast path (`release_gather`). During libMBD's
distributed matrix inverse, ScaLAPACK's `PDGETRF`/`PDGETRI` broadcast a column
panel described by a *strided vector* MPI datatype. On big-endian s390x, the
ch4 POSIX `release_gather` path miscomputes/serializes that datatype's byte
extent, so non-root ranks see a message size of `0` while the root advertises
`216`, and MPICH aborts the collective.

libMBD only issues standard ScaLAPACK calls here; it does not construct the MPI
datatype itself and cannot fix the broadcast. The defect belongs upstream in
MPICH. (Distros that build MPICH with the older `ch3` device — e.g. Debian and
Ubuntu — do **not** hit this, which is why it shows up only with `ch4` builds
such as Fedora's.)

## Reproduction

`devtools/reproduce-issue-76.sh` reproduces the crash in an emulated big-endian
s390x Fedora container (Fedora's MPICH is `ch4:ofi`). Requires Docker with QEMU
binfmt for s390x:

```sh
docker run --privileged --rm tonistiigi/binfmt --install s390x
devtools/reproduce-issue-76.sh          # crashes (expected, reproduces #76)
WORKAROUND=1 devtools/reproduce-issue-76.sh   # passes (workaround applied)
```

## Workaround (runtime MPICH setting)

Forcing MPICH's generic broadcast algorithm avoids the buggy ch4 POSIX
`release_gather` path. Set, at run time, before launching:

```sh
export MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM=mpir
```

With this set, `mbd_rsscs_deriv_expl` passes (numerical diff ~5e-9 vs the 1e-7
threshold). The broader hammer `MPIR_CVAR_DEVICE_COLLECTIVES=none` also works
but disables *all* device collectives (larger performance cost), so the
targeted CVAR above is preferred.

Settings that did **not** help (still crash): `MPIR_CVAR_ENABLE_INTRANODE_-
COLLECTIVES=0`, `MPIR_CVAR_BCAST_INTRANODE_ALGORITHM=mpir`,
`MPIR_CVAR_POSIX_BCAST_ALGORITHM=mpir`.

## Tested CVAR matrix (Fedora 40, MPICH ch4:ofi, s390x, 3 ranks)

| Setting                                          | Result |
| ------------------------------------------------ | ------ |
| _(baseline)_                                     | crash  |
| `MPIR_CVAR_DEVICE_COLLECTIVES=none`              | pass   |
| `MPIR_CVAR_ENABLE_INTRANODE_COLLECTIVES=0`       | crash  |
| `MPIR_CVAR_BCAST_INTRANODE_ALGORITHM=mpir`       | crash  |
| `MPIR_CVAR_POSIX_BCAST_ALGORITHM=mpir`           | crash  |
| `MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM=mpir`     | pass   |
