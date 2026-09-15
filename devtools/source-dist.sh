#!/bin/bash

set -ev
# GNU tar is `gtar` on macOS (where this is usually run by hand) and `tar` on
# Linux (e.g. in CI); both support --transform/--exclude, BSD tar does not. Pick
# the first GNU tar found, and fail with an actionable message rather than
# falling back to BSD tar and erroring later on an unknown --transform option.
TAR=$(command -v gtar || command -v tar)
if ! "${TAR}" --version 2>/dev/null | grep -q "GNU tar"; then
    echo "error: GNU tar is required (it provides --transform/--exclude)." >&2
    echo "Install it and re-run: 'brew install gnu-tar' on macOS, or" >&2
    echo "'apt-get install tar' on Debian/Ubuntu." >&2
    exit 1
fi
VERSION=$(git describe --tags --dirty=.dirty)
VERSION_FILE=cmake/libMBDVersionTag.cmake
SLUG=libmbd-${VERSION}
echo "set(VERSION_TAG ${VERSION})">${VERSION_FILE}
mkdir -p dist
ARCHIVE=dist/${SLUG}.tar.gz
"${TAR}" -vcz -f ${ARCHIVE} \
    --exclude "*pymbd*" --exclude "__pycache__" --exclude "conftest.py" --exclude ".*" \
    --transform "s,^,${SLUG}/," CMakeLists.txt cmake src tests LICENSE README.md
rm ${VERSION_FILE}
echo ${ARCHIVE}
