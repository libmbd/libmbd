#!/bin/bash

# Read the version and notes from CHANGELOG.md (Keep a Changelog format).
#
#   changelog.sh version         print the latest released version, e.g. 0.14.1
#   changelog.sh notes [VERSION] print the body of VERSION's section (default:
#                                the latest released version)

set -eu

CHANGELOG=${CHANGELOG:-CHANGELOG.md}

latest_version() {
    # First "## [X.Y.Z]" heading; "## [Unreleased]" has no digits so is skipped.
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
    # Lines between this version's heading and the next, blank edges trimmed.
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
