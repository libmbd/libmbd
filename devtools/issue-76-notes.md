Reference: <https://github.com/libmbd/libmbd/issues/76>

## Symptom (from the issue report)

Built with `-DENABLE_SCALAPACK_MPI=ON` and run under MPICH on s390x
(big-endian), the rSSCS derivative gradient tests abort. From #76
(`grad/mbd_rsscs_deriv_expl`, `BLACS grid: 1 x 3`, block size 3):

```
internal_Bcast(7723).......................: MPI_Bcast(... count=1, dtype=USER ...) failed
MPIDI_POSIX_mpi_bcast(238).................:
MPIDI_POSIX_mpi_release_gather_release(218): message sizes do not match
    across processes in the collective routine: Received 0 but expected 216
```

(`mbd_rsscs_deriv_impl_alpha` and `mbd_rsscs_deriv_impl_C6` fail identically.)
The same build passes on little-endian (x86-64, aarch64).

## Root cause — confirmed upstream, fixed in MPICH 5.0.1

This is a **big-endian bug in MPICH's release-gather collectives**, not a libMBD
bug. It was fixed upstream in **MPICH 5.0.1** (released 2026-04-10), whose
changelog reads:

> Fix bad cast in release-gather collectives that caused data loss issues on
> Big-Endian 64b arches (s390x)

That matches the failure on every axis: the `release_gather` path in the stack
trace, a bad cast / extent miscomputation, big-endian, s390x. The crash is the
`ch4` POSIX shared-memory broadcast (`MPIDI_POSIX_mpi_release_gather_release`)
mishandling the strided datatype that ScaLAPACK's `PDGETRF`/`PDGETRI` broadcast
during libMBD's distributed matrix inverse. libMBD only issues standard
ScaLAPACK calls here and cannot influence the broadcast.

Note the upstream changelog phrases the bug as *data loss*, whereas #76 shows a
hard size-mismatch abort. These are two manifestations of the same bad cast in
`release_gather` on big-endian — the corrupted extent can either truncate data
or trip the collective's size check.

## Fix

**Upgrade MPICH to ≥ 5.0.1.** That removes the bug at the source; no run-time
workaround is then needed.

## Workaround on older MPICH (< 5.0.1)

Force MPICH's generic broadcast algorithm so it bypasses the buggy `ch4` POSIX
`release_gather` path, at run time before launching:

```sh
export MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM=mpir
```

`MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM` is a real MPICH control variable
(declared in ch4's `shm/posix` collectives); selecting the `mpir` algorithm
routes the broadcast through the architecture-independent implementation instead
of the POSIX shared-memory fast path. The broader `MPIR_CVAR_DEVICE_COLLECTIVES=none`
also bypasses it but disables *all* device collectives, so the targeted CVAR is
preferred. Using a generic-algorithm CVAR to dodge a buggy device collective has
upstream precedent (e.g. pmodels/mpich #6983, #6984).

Only `ch4` builds with the POSIX shared-memory collectives have a
`release_gather` path, so only those are affected; the older `ch3` device has no
such path. `ch4` has been MPICH's default device since 3.4, so modern
distro builds (e.g. Fedora 40, `ch4:ofi`) take the affected path by default.

> Honest status: the specific run-time effect of the CVAR above on s390x has not
> been benchmarked here with reproducible logs — it rests on the CVAR's
> documented semantics and the now-confirmed upstream root cause, not on a
> measured pass/fail table. The decisive, verifiable fact is the MPICH 5.0.1 fix.

## Reproduction

`devtools/reproduce-issue-76.sh` builds and runs the failing gradient tests in an
emulated big-endian s390x Fedora container (Fedora's MPICH uses the `ch4`
device). It needs Docker with QEMU binfmt for s390x:

```sh
docker run --privileged --rm tonistiigi/binfmt --install s390x
devtools/reproduce-issue-76.sh                  # expected to reproduce on MPICH < 5.0.1
WORKAROUND=1 devtools/reproduce-issue-76.sh     # applies the CVAR workaround
```

The real reference environment remains native s390x (e.g. an IBM LinuxONE
instance or a distro s390x porterbox).

## References

- Issue #76 — <https://github.com/libmbd/libmbd/issues/76>
- MPICH 5.0.1 release note — <https://www.mpich.org/2026/04/10/mpich-5-0-1-released/>
- MPICH releases / changelog — <https://github.com/pmodels/mpich/releases>
- MPICH `CHANGES` (ch4 default since 3.4) — <https://github.com/pmodels/mpich/blob/main/CHANGES>
- CVAR-as-workaround precedent — <https://github.com/pmodels/mpich/issues/6983>, <https://github.com/pmodels/mpich/issues/6984>
