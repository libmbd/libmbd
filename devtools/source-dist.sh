#!/bin/bash

set -ev
VERSION=$(git describe --tags --dirty=.dirty)
VERSION_FILE=cmake/libMBDVersionTag.cmake
SLUG=libmbd-${VERSION}
echo "set(VERSION_TAG ${VERSION})">${VERSION_FILE}
mkdir -p dist
ARCHIVE=dist/${SLUG}.tar.gz
gtar -vcz -f ${ARCHIVE} \
    --exclude "*pymbd*" --exclude "__pycache__" --exclude "conftest.py" --exclude ".*" \
    --transform "s,^,${SLUG}/," CMakeLists.txt cmake src tests LICENSE README.md
rm ${VERSION_FILE}
echo ${ARCHIVE}
