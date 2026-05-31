# Releasing

The version is derived automatically from the latest Git tag (via
poetry-dynamic-versioning and CMake), so there is no version string to bump in
the source tree. A release is a changelog update followed by a tag.

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
2. Merge the release PR.
3. Tag the merge commit and push the tag:

   ```
   git tag X.Y.Z
   git push origin X.Y.Z
   ```

   The tag drives the package version automatically.
