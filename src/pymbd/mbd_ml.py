# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

from pymbd import from_volumes
from pymbd.fortran import MBDGeom

try:
    import ase.io
    from ase.units import Bohr
except ImportError as e:
    raise ImportError(
        'pymbd.mbd_ml requires ase, which is not installed. See '
        'requirements-mbd-ml.txt.'
    ) from e

#: How to get the optional dependencies that are not on PyPI.
_INSTALL_HINT = (
    'pymbd.mbd_ml requires so3lr, which is not on PyPI. Install it with\n\n'
    '    pip install -e . -r requirements-mbd-ml.txt\n\n'
    'from a checkout, or with\n\n'
    '    pip install pymbd -r '
    'https://raw.githubusercontent.com/libmbd/libmbd/master/requirements-mbd-ml.txt'
)


def _so3lr():
    """Import so3lr on demand, with a message pointing at how to install it.

    so3lr is imported lazily rather than at module level so that the parts of
    this module that do not need it -- the periodicity checks and the stress
    conversion -- stay importable, and testable, without it.
    """
    try:
        from so3lr.cli.so3lr_eval import evaluate_so3lr_on
        from so3lr.models import MBD_ML_MODELS, model_path
    except ImportError as e:
        raise ImportError(_INSTALL_HINT) from e
    return evaluate_so3lr_on, MBD_ML_MODELS, model_path


def _shipped_models():
    """Names of the MBD-ML models so3lr ships, or () if so3lr is unavailable."""
    try:
        return _so3lr()[1]
    except ImportError:
        return ()


DEFAULT_MODEL = 'sv2j_b64_l2d_42e_16hh_10_24novv'


def resolve_model(model):
    """Return the directory holding the MBD-ML model `model`.

    `model` is either the name of a model shipped by so3lr (see
    `so3lr.models.MBD_ML_MODELS`) or the path of a directory holding one.
    """
    shipped = _shipped_models()
    if model in shipped:
        return _so3lr()[2](model)
    path = Path(model).expanduser().resolve()
    if (path / 'hyperparameters.json').is_file():
        return path
    if not shipped:
        # so3lr is missing, so `model` could well be one of its model names
        # rather than a bad path; the install hint is the more useful error.
        _so3lr()
    raise FileNotFoundError(
        f'{model!r} is neither an MBD-ML model shipped by so3lr '
        f'({", ".join(shipped)}) nor a directory holding one'
    )


def ratios_from_mbdml(atoms, model=DEFAULT_MODEL):
    """Predict the alpha_0 and C6 ratios of a structure with the MBD-ML model.

    :param atoms: ASE ``Atoms`` object
    :param model: name of an MBD-ML model shipped by so3lr, or the path of a
        directory holding one

    Returns a dict with the per-atom ``'a0'`` and ``'c6'`` ratios.
    """
    model_path = resolve_model(model)
    tmpdir = tempfile.mkdtemp(prefix='mbdml-')
    mbdml_in_filename = os.path.join(tmpdir, 'mbdml_in.extxyz')
    mbdml_out_filename = os.path.join(tmpdir, 'mbdml_out.extxyz')

    try:
        ase.io.write(
            mbdml_in_filename,
            atoms,
            format='extxyz',
            write_info=True,
            write_results=True,
        )

        evaluate_so3lr_on = _so3lr()[0]
        _ = evaluate_so3lr_on(
            datafile=mbdml_in_filename,
            batch_size=1,
            lr_cutoff=0.1,
            dispersion_damping=2.0,
            jit_compile=False,
            save_to=mbdml_out_filename,
            model_path=model_path,
            precision='float32',
            targets='hirshfeld_ratios,c6_ratios',
            log_file=None,
        )

        atoms_eval = ase.io.read(mbdml_out_filename, format='extxyz')
        c6 = atoms_eval.arrays['c6_ratios_so3lr']
        a0 = atoms_eval.arrays['hirshfeld_ratios_so3lr']
    except Exception:
        print(f'MBD-ML evaluation failed, temporary files kept in {tmpdir}')
        raise

    # Remove temporary xyz files, as otherwise so3lr eval fails in the second step
    shutil.rmtree(tmpdir)

    combined_ratios = np.concatenate([c6, a0])
    ratio_min = 0.05
    ratio_max = 3.0
    if np.any((combined_ratios < ratio_min) | (combined_ratios > ratio_max)):
        print(
            f"\n{'!' * 50}\nWARNING: a0 or c6 ratios outside [{ratio_min}, {ratio_max}]!\nThis indicates that either your system is pathological or that the MBD-ML is\nextrapolating and th\
at the result is potentially not reliable. Proceed with care!\n{'!' * 50}"
        )

    return {'c6': c6, 'a0': a0}


def compute_stress_from_lattice_gradient(
    lattice, coords_cartesian, dE_dlattice, dE_dcoords
):
    """Convert a lattice gradient into a stress tensor (a.u.).

    :param lattice: lattice vectors as rows
    :param coords_cartesian: Cartesian atomic coordinates
    :param dE_dlattice: energy gradient with respect to the lattice vectors
    :param dE_dcoords: energy gradient with respect to the coordinates
    """
    term1 = lattice.T @ dE_dlattice
    term2 = coords_cartesian.T @ dE_dcoords

    stress_times_volume = term1 + term2
    cell_vol = abs(np.linalg.det(lattice))

    return stress_times_volume / cell_vol


def mbd_properties_from_structure(atoms, beta, k_grid=None, model=DEFAULT_MODEL):
    """Compute the MBD energy, forces and stress of a structure with MBD-ML ratios.

    :param atoms: ASE ``Atoms`` object
    :param float beta: MBD range-separation parameter
    :param k_grid: k-point grid, required for periodic systems
    :param model: name of an MBD-ML model shipped by so3lr, or the path of a
        directory holding one

    Returns a dict with the energy ``'E'``, forces ``'F'`` and, for periodic
    systems, the stress ``'S'``, in atomic units.
    """

    if any(atoms.pbc) and k_grid is None:
        raise ValueError('k_grid must be given for periodic systems')

    ratios_dict = ratios_from_mbdml(atoms, model=model)

    a0_ratios = ratios_dict['a0']
    C6_ratios = ratios_dict['c6']

    atom_pos = atoms.get_positions() / Bohr
    atom_species = atoms.get_chemical_symbols()
    n_atoms = len(atoms)

    if any(atoms.pbc):
        lattice_vecs = atoms.cell[:, :] / Bohr
    else:
        lattice_vecs = None

    a0_free, C6_free, _ = from_volumes(atom_species, np.ones(n_atoms))
    Rvdw = 2.5 * a0_free ** (1 / 7) * a0_ratios ** (1 / 3)

    a0 = a0_free * a0_ratios
    C6 = C6_free * C6_ratios

    confMBD = MBDGeom(coords=atom_pos, lattice=lattice_vecs, k_grid=k_grid)

    # mbd_energy() returns energies in Ha and energy gradients in Ha/Bohr
    if any(atoms.pbc):
        E, gradE, dE_dL = confMBD.mbd_energy(
            a0,
            C6,
            R_vdw=Rvdw,
            beta=beta,
            damping='fermi,dip',
            variant='plain',
            force=True,
        )
    else:
        E, gradE = confMBD.mbd_energy(
            a0,
            C6,
            R_vdw=Rvdw,
            beta=beta,
            damping='fermi,dip',
            variant='plain',
            force=True,
        )

    F = (-1.0) * gradE

    if any(atoms.pbc):
        S = compute_stress_from_lattice_gradient(lattice_vecs, atom_pos, dE_dL, -F)
        return {'E': E, 'F': F, 'S': S}
    else:
        return {'E': E, 'F': F}
