# Releasing

The version is derived automatically from the latest Git tag (via
poetry-dynamic-versioning and CMake), so there is no version string to bump in
the source tree. A release is driven entirely by the changelog: the topmost
`## [X.Y.Z]` heading in `CHANGELOG.md` names the version, and merging a PR
labelled `release` tells CI to tag that commit and publish everything.

The changelog (`CHANGELOG.md`) is maintained by hand following
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/). Add notable,
user-facing changes under the `## [Unreleased]` heading as part of the PR that
makes them.

To cut release `X.Y.Z`:

1. Open a release PR that edits `CHANGELOG.md`:
   - Rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD`.
   - Add a fresh, empty `## [Unreleased]` heading above it.
   - Update the link references at the bottom: point `[unreleased]` at
     `compare/X.Y.Z...HEAD` and add `[X.Y.Z]: .../compare/<prev>...X.Y.Z`.
2. Add the **`release`** label to the PR.
3. Merge it. That is the entire trigger — CI does the rest.

No tag is created by hand: the `Release` workflow reads `X.Y.Z` from the
changelog and creates the tag itself, on the exact merge commit.

## What CI does on release

Merging a `release`-labelled PR runs the `Release` workflow
(`.github/workflows/release.yaml`), which:

1. **prepare** — reads the version from the top of `CHANGELOG.md`
   (`devtools/changelog.sh version`), refuses to proceed if that tag already
   exists, and creates and pushes the tag on the merge commit. The release
   notes are extracted from that changelog section
   (`devtools/changelog.sh notes`).
2. **build** — reuses the `Build distributions` workflow
   (`.github/workflows/build-dist.yaml`, the same one that validates packaging
   on pull requests) against the new tag, building the pyMBD source
   distribution and the libMBD source archive.
3. **publish** — creates the GitHub Release **as a draft**, attaches both source
   archives while it is still mutable, uploads the pyMBD sdist to
   [PyPI](https://pypi.org/project/pymbd/) via
   [trusted publishing](https://docs.pypi.org/trusted-publishers/) (no API token
   is stored; the `pymbd` project trusts this repository's `release.yaml`
   workflow under the `pypi` environment), and only then flips the release to
   published.

Building the assets and attaching them to a draft *before* publishing is what
makes the flow compatible with
[immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases):
an immutable release is frozen the moment it is published, so every asset must
already be in place. The previous flow published the Release first and attached
assets afterwards, which immutability forbids.

Only the source distribution goes to PyPI: pyMBD's C extension is compiled
against an external libMBD at install time, so a binary wheel would bake in one
libMBD build and break for everyone else. The libMBD source archive
(`libmbd-X.Y.Z.tar.gz`, built by `devtools/source-dist.sh`) is published only as
a Release asset, alongside GitHub's auto-generated source code archives.

Conda-forge's autotick bot then picks up the new release and opens a
version-bump PR on the
[`libmbd-feedstock`](https://github.com/conda-forge/libmbd-feedstock), which
builds and publishes the `libmbd` conda package once merged.

## One-time setup

Two repository settings need to be aligned with this workflow:

- **PyPI trusted publisher.** The `pymbd` project on PyPI must trust the
  workflow file **`release.yaml`** (environment `pypi`). If it still lists the
  old `publish.yaml`, update it under the project's Publishing settings, or the
  OIDC upload will be rejected.
- **Immutable releases.** Enabling immutable releases (repository → Settings →
  General, or at the organization level) is optional but is the whole point of
  the draft-then-publish flow above; the workflow already produces releases that
  satisfy it.
