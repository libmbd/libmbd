# Releasing

The version is derived automatically from the latest Git tag (via
poetry-dynamic-versioning and CMake), so there is no version string to bump in
the source tree. A release is a changelog update followed by a tagged GitHub
Release, which CI turns into the published packages and assets.

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
2. Merge the release PR (its CI run also validates the release commit).
3. Publish a [GitHub Release](https://github.com/libmbd/libmbd/releases/new)
   targeting the merge commit, using `X.Y.Z` as the tag. GitHub creates the tag
   as part of publishing the Release, so no local `git`/`git push` is needed.
   The tag must be the bare version number (matching the
   `poetry-dynamic-versioning` pattern); it is what drives the package version.

## What CI does on release

Publishing the Release triggers the `Build distributions` workflow
(`.github/workflows/build-dist.yaml`) — the same workflow that validates
packaging changes on pull requests — which builds the pyMBD source distribution
and the libMBD source archive and stores them as run artifacts. When that build
succeeds, the `Publish` workflow (`.github/workflows/publish.yaml`) picks it up
via `workflow_run` (only acting on a successful build that came from a published
Release) and:

- uploads the pyMBD source distribution to
  [PyPI](https://pypi.org/project/pymbd/) via
  [trusted publishing](https://docs.pypi.org/trusted-publishers/) — no API token
  is stored; the `pymbd` project trusts this repository's `publish.yaml` workflow
  under the `pypi` environment; and
- attaches both source archives to the Release — the pyMBD sdist
  (`pymbd-X.Y.Z.tar.gz`) and the libMBD source archive (`libmbd-X.Y.Z.tar.gz`,
  built by `devtools/source-dist.sh`) — alongside GitHub's auto-generated source
  code archives. These previously had to be built and uploaded by hand.

Only the source distribution goes to PyPI: pyMBD's C extension is compiled
against an external libMBD at install time, so a binary wheel would bake in one
libMBD build and break for everyone else. The libMBD source archive is published
only as a Release asset.

Conda-forge's autotick bot then picks up the new release and opens a
version-bump PR on the
[`libmbd-feedstock`](https://github.com/conda-forge/libmbd-feedstock), which
builds and publishes the `libmbd` conda package once merged.
