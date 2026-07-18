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
3. Merge it. That builds and stages a **draft** GitHub Release.
4. Review the draft Release, then **publish it by hand**. Publishing is what
   uploads the package to PyPI.

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
3. **draft** — creates the GitHub Release **as a draft** with both source
   archives attached.

The workflow stops there. A maintainer reviews the draft — the changelog notes,
the attached archives — and publishes it by hand. Publishing fires the
`release: published` event, which runs the `Publish` workflow
(`.github/workflows/publish.yaml`): it downloads the pyMBD sdist back from the
Release and uploads it to [PyPI](https://pypi.org/project/pymbd/) via
[trusted publishing](https://docs.pypi.org/trusted-publishers/) (no API token is
stored; the `pymbd` project trusts this repository's `publish.yaml` workflow
under the `pypi` environment). PyPI therefore receives the exact artifact
attached to the Release.

Attaching every asset to the draft *before* it is published is what makes the
flow compatible with
[immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases):
an immutable release is frozen the moment it is published, so every asset must
already be in place. Publishing to PyPI reads the Release but does not modify it.
The previous flow published the Release first and attached assets afterwards,
which immutability forbids.

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

- **Immutable releases.** Enabling immutable releases (repository → Settings →
  General, or at the organization level) is optional but is the whole point of
  the draft-then-publish flow above; the workflow already produces releases that
  satisfy it.

The PyPI trusted publisher is unchanged: the `pymbd` project still trusts the
`publish.yaml` workflow under the `pypi` environment, exactly as before.
