import os
import datetime
from configparser import ConfigParser

conf = ConfigParser()
conf.read('../setup.cfg')

project = 'libmbd'
version = conf.get('libmbd:info', 'version')
author = conf.get('libmbd:info', 'author')
description = conf.get('libmbd:info', 'description')

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'breathe',
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
    'github_user': 'azag0',
    'github_repo': 'libmbd',
    'badge_branch': 'master',
    'codecov_button': True,
    'travis_button': True,
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


def setup(app):
    app.connect('autodoc-skip-member', skip_namedtuples)
