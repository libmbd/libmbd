# Notes for Claude

## Changelog

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and is maintained by hand. Add an entry for a user-visible change as part of the
change itself, not as a follow-up:

- Put it under `## [Unreleased]`, in one of the categories the file already uses
  — `Added`, `Changed`, `Removed`, `Fixed` — adding the `###` heading if that
  category is not there yet.
- Link the pull request or issue it came from:
  `- Support for cffi 2.x ([#126](https://github.com/libmbd/libmbd/pull/126))`
- Say what changed from the outside, for someone deciding whether an upgrade
  affects them, rather than how it was implemented. Prefix a breaking change
  with `**Breaking:**`.

**Never add a `## [X.Y.Z]` version heading.** The release workflow reads the
version from the topmost one and tags a release from it, so inventing a heading
starts a release. Turning `[Unreleased]` into a version is the release process's
job — see `RELEASING.md`.

Changes with no effect a user could observe — refactoring, tests, CI,
documentation, build plumbing — get no entry. Leave them out rather than padding
the list.
