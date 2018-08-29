#!/usr/bin/env python3
import datetime
import subprocess
import os


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'breathe',
]
source_suffix = '.rst'
master_doc = 'index'
project = 'libmbd'
author = 'Jan Hermann'
copyright = f'2018-{datetime.date.today().year}, {author}'
version = '0.4.0a1'
release = version
language = None
exclude_patterns = ['build', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = True
html_theme = 'alabaster'
html_theme_options = {
    'description': 'Many-body dispersion library',
    'github_button': True,
    'github_user': 'azag0',
    'github_repo': 'libmbd',
    'extra_nav_links': {
        'Fortran API': 'doxygen/namespacembd.html',
        'Core Fortran module': 'doxygen/namespacembd__core.html',
    },
}
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
    ]
}
htmlhelp_basename = f'{project}doc'

breathe_projects_source = {
    'fortran-c-api': ('../src', ['mbd.h'])
}


def skip_namedtuples(app, what, name, obj, skip, options):
    if hasattr(obj, '_source'):
        return True


def run_doxygen(app, exception):
    if os.environ.get('READTHEDOCS') != 'True':
        return
    if os.path.exists('_build/html/oxygen'):
        return
    subprocess.check_call(
        '(cat Doxyfile; echo "OUTPUT_DIRECTORY=_build/html") | doxygen -',
        shell=True
    )


def setup(app):
    app.connect('autodoc-skip-member', skip_namedtuples)
    app.connect('build-finished', run_doxygen)
