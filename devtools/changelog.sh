#!/bin/bash

# Read release metadata from CHANGELOG.md (Keep a Changelog format).
#
# The changelog is the single source of truth for what a release *is*: the
# topmost versioned section (skipping the [Unreleased] heading) names the
# version to tag, and its body becomes the GitHub Release notes. This lets a
# changelog-only PR drive the whole release, with the version living in a
# reviewed file rather than being typed into a tag by hand.
#
# Usage:
#   changelog.sh version         print the latest released version, e.g. 0.14.1
#   changelog.sh notes [VERSION] print the body of VERSION's section, defaulting
#                                to the latest released version

set -eu

CHANGELOG=${CHANGELOG:-CHANGELOG.md}

latest_version() {
    # First heading of the form "## [X.Y.Z]". "## [Unreleased]" carries no
    # digits, so the version pattern skips it and lands on the newest release.
    sed -nE 's/^## \[([0-9]+\.[0-9]+\.[0-9]+)\].*/\1/p' "$CHANGELOG" | head -n1
}

case "${1:-}" in
version)
    version=$(latest_version)
    if [ -z "$version" ]; then
        echo "error: no released version found in $CHANGELOG" >&2
        exit 1
    fi
    echo "$version"
    ;;
notes)
    version=${2:-$(latest_version)}
    if [ -z "$version" ]; then
        echo "error: no released version found in $CHANGELOG" >&2
        exit 1
    fi
    # Everything between this version's heading and the next "## [" heading,
    # with leading and trailing blank lines trimmed so the release notes start
    # and end on real content.
    awk -v ver="$version" '
        index($0, "## [" ver "]") == 1 { inside = 1; next }
        inside && /^## \[/            { exit }
        inside                        { buf[++n] = $0 }
        END {
            start = 1; while (start <= n && buf[start] ~ /^[[:space:]]*$/) start++
            end = n; while (end >= start && buf[end] ~ /^[[:space:]]*$/) end--
            for (i = start; i <= end; i++) print buf[i]
        }
    ' "$CHANGELOG"
    ;;
*)
    echo "usage: changelog.sh {version|notes [VERSION]}" >&2
    exit 1
    ;;
esac
