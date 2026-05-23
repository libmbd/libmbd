#!/bin/bash
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
    exit 0
fi

# Install system dependencies for libMBD (Fortran compiler + LAPACK/BLAS)
sudo apt-get install -yq --no-install-suggests --no-install-recommends \
    gfortran libblas-dev liblapack-dev

# Set up Python virtual environment
VENV="$HOME/env"
if [ ! -d "$VENV" ]; then
    python3 -m venv "$VENV"
fi
source "$VENV/bin/activate"
pip install -U pip --quiet

# Build and install libMBD (Fortran) + pyMBD (Python bindings + test deps)
cd "$CLAUDE_PROJECT_DIR"
export VIRTUAL_ENV="$VENV"
export LIBMBD_PREFIX="$VENV"
make install_editable

# Persist environment variables for the rest of the session
echo "export VIRTUAL_ENV=$VENV" >> "$CLAUDE_ENV_FILE"
echo "export PATH=$VENV/bin:\$PATH" >> "$CLAUDE_ENV_FILE"
echo "export LIBMBD_PREFIX=$VENV" >> "$CLAUDE_ENV_FILE"
