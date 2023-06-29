import datetime
import os
import subprocess
import sys

import toml

sys.path.insert(0, os.path.abspath('../src'))


with open('../pyproject.toml') as f:
    metadata = toml.load(f)['tool']['poetry']

project = 'libMBD'
release = version = (
    subprocess.run(['poetry', 'version'], capture_output=True, cwd='..')
    .stdout.decode()
    .split()[1]
)
author = ' '.join(metadata['authors'][0].split()[:-1])
description = metadata['description']

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}
source_suffix = '.rst'
master_doc = 'index'
copyright = f'2018-{datetime.date.today().year}, {author}'
language = None
exclude_patterns = ['build', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = True
html_theme = 'alabaster'
html_theme_options = {
    'description': description,
    'github_button': True,
    'github_user': 'libmbd',
    'github_repo': 'libmbd',
    'badge_branch': 'master',
    'codecov_button': True,
}
html_sidebars = {
    '**': ['about.html', 'navigation.html', 'relations.html', 'searchbox.html']
}

autodoc_default_options = {'members': True}
autodoc_inherit_docstrings = False
autodoc_mock_imports = [
    'numpy',
    'scipy',
    'pymbd._libmbd',
    'mpi4py',
]


def skip_namedtuples(app, what, name, obj, skip, options):
    if hasattr(obj, '_source'):
        return True


def setup(app):
    app.connect('autodoc-skip-member', skip_namedtuples)
