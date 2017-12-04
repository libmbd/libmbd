#!/usr/bin/env python -B
from __future__ import print_function
import os
import json
import sys
from argparse import ArgumentParser

conf_file = 'conf.json'
sources = [
    'mbd.f90',
    'mbd_c_api.f90',
    'mbd_helper.f90',
    'mbd_helper_dev.f90',
    'mbd_interface.f90',
    'mbd_math.f90',
    'mbd_repulsion.f90',
]


def configure(opts, env):
    blddir = opts.blddir
    srcdir = os.path.realpath(opts.srcdir)
    env.update(dict(
        VPATH=srcdir,
        FCOMPILEOPTS='',
        FOBJS=[fname.replace('.f90', '.o') for fname in sources],
    ))
    fortran_tasks = {fname: {
        'source': os.path.join(srcdir, fname),
        'args': env['FC'] + ['-c', '-o', fname.replace('.f90', '.o')] + env['FFLAGS'],
    } for fname in sources}
    if not os.path.isdir(blddir):
        os.mkdir(blddir)
    with open('Makefile.build') as f:
        makefile = f.read()
    with open(os.path.join(blddir, 'Makefile'), 'w') as f:
        for varname, value in env.items():
            f.write('{0} = {1}\n'.format(
                varname,
                ' '.join(value) if isinstance(value, list) else str(value)
            ))
        f.write('\n')
        f.write(makefile)
    with open(os.path.join(blddir, 'fcompile.json'), 'w') as f:
        json.dump(fortran_tasks, f)
    with open('Makefile.in') as f:
        makefile = f.read()
    with open('Makefile', 'w') as f:
        f.write('BLDDIR = {0}\n\n'.format(blddir))
        f.write(makefile)


def save_conf(opts, argv):
    quoted_argv = ["'{0}'".format(arg) if ' ' in arg else arg for arg in argv]
    with open(conf_file, 'w') as f:
        json.dump({
            'argline': ' '.join(quoted_argv),
            'argv': argv,
            'blddir': opts.blddir,
        }, f, indent=4)


def main(argv=sys.argv[1:]):
    parser = ArgumentParser(
        # usage='usage: ./configure.py [options] FC=<fortran compiler> [VAR=VAL...]'
    )
    arg = parser.add_argument
    arg('--no-fcompile', action='store_true', help='options for fcompile')
    arg('--srcdir', default='src', help='source directory')
    arg('--blddir', default='_build', help='object file directory')
    arg('--prefix', default='/usr/local', help='prefix path')
    arg('VAR', nargs='*', help='shared library directory')
    if not argv and os.path.exists(conf_file):
        with open(conf_file) as f:
            argv = json.load(f)['argv']
    args = parser.parse_args(argv)
    env = {var: val.split() for var, val in (arg.split('=', 1) for arg in args.VAR)}
    save_conf(args, argv)
    print('Configuring with: {0!r}'.format(argv))
    configure(args, env)


if __name__ == '__main__':
    main()
