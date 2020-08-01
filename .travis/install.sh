if [[ $TRAVIS_OS_NAME = osx ]]; then
    python3 -m venv ~/virtualenv
    source ~/virtualenv/bin/activate
    python --version
    pip --version
fi
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
source $HOME/.poetry/env
pip install codecov pytest-cov "coverage<5"
echo VIRTUAL_ENV=$VIRTUAL_ENV
