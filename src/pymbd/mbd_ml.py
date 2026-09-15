# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import shutil
import tempfile
from pymbd import from_volumes
from pymbd.fortran  import MBDGeom
import numpy as np

try:
    import ase.io
    from ase.units import Bohr, Hartree
except ImportError:
    print('ase package cannot be imported but is required to use MBD-ML. Please install.')
    sys.exit(1)

try:
    from so3lr.cli.so3lr_eval import evaluate_so3lr_on
except ImportError:
    print('so3lr package cannot be imported. Please reinstall pymbd via pip install ".[mbd-ml]"')
    sys.exit(1)

def ratios_from_mbdml(atoms):
    '''This function is a wrapper for the MBD-ML model, which predicts a0 and C6 ratios given a molecular or crystal structure'''
    tmpdir = tempfile.mkdtemp(prefix='mbdml-')
    mbdml_in_filename = os.path.join(tmpdir, 'mbdml_in.extxyz')
    mbdml_out_filename = os.path.join(tmpdir, 'mbdml_out.extxyz')

    try:
        ase.io.write(mbdml_in_filename, atoms, format='extxyz', write_info=True, write_results=True)

        #model_path=f"{os.path.dirname(os.path.realpath(__file__))}/mbd_ml_model/sv2j_b128_l2d_42e_16hh_10_24novv"
        model_path=f"{os.path.dirname(os.path.realpath(__file__))}/mbd_ml_model/sv2j_b64_l2d_42e_16hh_10_24novv"
        print(f'Path of model: {model_path}')
        _ = evaluate_so3lr_on(
                datafile = mbdml_in_filename,
                batch_size = 1,
                lr_cutoff = 0.1,
                dispersion_damping = 2.0,
                jit_compile = False,
                save_to = mbdml_out_filename,
                model_path = model_path,
                precision = "float32",
                targets = "hirshfeld_ratios,c6_ratios",
                log_file = None
                )

        atoms_eval = ase.io.read(mbdml_out_filename, format='extxyz')
        c6 = atoms_eval.arrays["c6_ratios_so3lr"]
        a0 = atoms_eval.arrays["hirshfeld_ratios_so3lr"]
    except Exception:
        print(f'MBD-ML evaluation failed, temporary files kept in {tmpdir}')
        raise

    #Remove temporary xyz files, as otherwise so3lr eval fails in the second step
    shutil.rmtree(tmpdir)

    combined_ratios = np.concatenate([c6, a0])
    ratio_min = 0.05
    ratio_max = 3.0
    if np.any((combined_ratios < ratio_min) | (combined_ratios > ratio_max)):
        print(f"\n{'!'*50}\nWARNING: a0 or c6 ratios outside [{ratio_min}, {ratio_max}]!\nThis indicates that either your system is pathological or that the MBD-ML is\nextrapolating and th\
at the result is potentially not reliable. Proceed with care!\n{'!'*50}")

    return {'c6' : c6, 'a0': a0}

def compute_stress_from_lattice_gradient(lattice, coords_cartesian, dE_dlattice, dE_dcoords):
    term1 = lattice.T @ dE_dlattice
    term2 = coords_cartesian.T @ dE_dcoords

    stress_times_volume = term1 + term2
    cell_vol = abs(np.linalg.det(lattice))

    return stress_times_volume/cell_vol



def mbd_properties_from_structure(atoms, beta, k_grid=None):
    '''Given a molecular or crystal structure, this function uses ratios_from_mbdml to obtain ratios and then computes energy, forces and stress'''

    if any(atoms.pbc) and k_grid is None:
        raise ValueError(
                "k_grid must be given for periodic systems"
        )

    ratios_dict = ratios_from_mbdml(atoms)

    a0_ratios = ratios_dict['a0']
    C6_ratios = ratios_dict['c6']

    atom_pos = atoms.get_positions() / Bohr
    atom_species = atoms.get_chemical_symbols()
    atom_Z = atoms.get_atomic_numbers()
    n_atoms = len(atoms)

    if any(atoms.pbc):
        lattice_vecs = atoms.cell[:,:] / Bohr
    else:
        lattice_vecs = None

    a0_free, C6_free, _ = from_volumes(atom_species, np.ones(n_atoms))
    Rvdw = 2.5 * a0_free**(1/7) * a0_ratios**(1/3)

    a0 = a0_free * a0_ratios
    C6 = C6_free * C6_ratios

    confMBD = MBDGeom(coords=atom_pos, lattice=lattice_vecs, k_grid=k_grid)

    #mbd_energy() returns energies in Ha and energy gradients in Ha/Bohr
    if any(atoms.pbc):
        E,gradE,dE_dL = confMBD.mbd_energy(a0,C6,R_vdw=Rvdw,beta=beta,damping='fermi,dip',variant='plain', force=True)
    else:
        E,gradE = confMBD.mbd_energy(a0,C6,R_vdw=Rvdw,beta=beta,damping='fermi,dip',variant='plain', force=True)

    F = (-1.0) * gradE

    if any(atoms.pbc):
        S = compute_stress_from_lattice_gradient(lattice_vecs, atom_pos, dE_dL, -F)
        return {'E' : E, 'F' : F, 'S' : S}
    else:
        return {'E' : E, 'F' : F}


