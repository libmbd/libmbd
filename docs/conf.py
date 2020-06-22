import datetime
import os
import sys
from unittest.mock import MagicMock

import toml

sys.path.insert(0, os.path.abspath('../src'))


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()


MOCK_MODULES = [
    'numpy',
    'numpy.linalg',
    'numpy.polynomial',
    'numpy.polynomial.legendre',
    'scipy',
    'scipy.special',
    'pymbd._libmbd',
    'mpi4py',
]
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

with open('../pyproject.toml') as f:
    metadata = toml.load(f)['tool']['poetry']

project = 'Libmbd'
version = metadata['version']
author = ' '.join(metadata['authors'][0].split()[:-1])
description = metadata['description']

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]
source_suffix = '.rst'
master_doc = 'index'
copyright = f'2018-{datetime.date.today().year}, {author}'
release = version
language = None
exclude_patterns = ['build', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = True
html_theme = 'alabaster'
html_theme_options = {
    'description': description,
    'github_button': True,
    'github_user': 'jhrmnn',
    'github_repo': 'libmbd',
    'badge_branch': 'master',
    'codecov_button': True,
    'travis_button': True,
}
html_sidebars = {
    '**': ['about.html', 'navigation.html', 'relations.html', 'searchbox.html']
}
htmlhelp_basename = f'{project}doc'


def skip_namedtuples(app, what, name, obj, skip, options):
    if hasattr(obj, '_source'):
        return True


def setup(app):
    app.connect('autodoc-skip-member', skip_namedtuples)
