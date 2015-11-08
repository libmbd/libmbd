#!/usr/bin/env python3
import json
from mpi4py import MPI
import sys
from pathlib import Path


bohr = 0.529177249
ntasks = MPI.COMM_WORLD.Get_size()
myid = MPI.COMM_WORLD.Get_rank()

free_atom_db = json.load((Path(__file__).parent/'free_atoms.json').open())


def get_free_atom_data(species):
    return list(zip(*[(at['alpha_0'], at['C6'], at['R_vdw']) for at in
                      [free_atom_db[sp] for sp in species]]))


def printout(s, each=False):
    if each or myid == 0:
        sys.stdout.write('{}\n'.format(s))


def printerr(s, each=False):
    if each or myid == 0:
        sys.stderr.write('{}\n'.format(s))


class ArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        try:
            return obj.tolist()
        except AttributeError:
            return super().default(obj)
