#!/bin/bash

set -ev
# GNU tar is `gtar` on macOS (where this is usually run by hand) and `tar` on
# Linux (e.g. in CI); both support --transform/--exclude, BSD tar does not.
TAR=$(command -v gtar || command -v tar)
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
