# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json
import os
import shutil
import tempfile
import types
import warnings
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
    """Import so3lr and its jax stack on demand, with an install hint.

    so3lr is imported lazily rather than at module level so that the parts of
    this module that do not need it -- the periodicity checks and the stress
    conversion -- stay importable, and testable, without it.
    """
    try:
        import jax
        import jraph
        from mlff.data import AseDataLoaderSparse
        from mlff.mdx.potential.mlff_potential_sparse import load_model_from_workdir
        from mlff.utils import jraph_utils
        from so3lr.models import MBD_ML_MODELS, model_path
    except ImportError as e:
        raise ImportError(_INSTALL_HINT) from e
    return types.SimpleNamespace(
        jax=jax,
        jraph=jraph,
        AseDataLoaderSparse=AseDataLoaderSparse,
        load_model_from_workdir=load_model_from_workdir,
        jraph_utils=jraph_utils,
        MBD_ML_MODELS=MBD_ML_MODELS,
        model_path=model_path,
    )


#: Long-range settings the MBD-ML models were evaluated with. They have no
#: long-range terms of their own (both ``*_energy_bool`` are false), but the
#: loader insists on values.
LR_CUTOFF, LR_DAMPING = 0.1, 2.0

_MODELS = {}


def _load_model(model):
    """Load a model once and cache it, keyed by its resolved directory.

    Loading the checkpoint dominates a single evaluation, so a module-level
    cache is what makes repeated calls -- a relaxation, a benchmark -- cheap.
    """
    path = str(resolve_model(model))
    if path not in _MODELS:
        so3lr = _so3lr()
        with open(os.path.join(path, 'hyperparameters.json')) as f:
            cutoff = json.load(f)['model']['cutoff']
        net, params = so3lr.load_model_from_workdir(
            path,
            model='so3krates',
            from_file=False,
            long_range_kwargs={
                'cutoff_lr': LR_CUTOFF,
                'dispersion_energy_cutoff_lr_damping': LR_DAMPING,
                'neighborlist_format_lr': 'sparse',
            },
        )
        _MODELS[path] = (net, params, cutoff)
    return _MODELS[path]


def _graph_inputs(atoms, cutoff):
    """Build the padded graph the model expects from an ASE ``Atoms``.

    The neighbour lists come from mlff's own loader, which reads a file, so the
    structure goes through a temporary extxyz. Everything after this is
    in-process, which is what makes the ratios differentiable.
    """
    so3lr = _so3lr()
    tmpdir = tempfile.mkdtemp(prefix='mbdml-')
    try:
        path = os.path.join(tmpdir, 'mbdml_in.extxyz')
        ase.io.write(path, atoms, format='extxyz', write_info=True, write_results=True)
        data, _ = so3lr.AseDataLoaderSparse(path).load(
            cutoff=cutoff,
            cutoff_lr=LR_CUTOFF,
            calculate_neighbors_lr=False,
            pick_idx=np.arange(1),
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    batch = so3lr.jraph.batch_np([data[0]])
    n_pairs = int(getattr(batch, 'n_pairs', np.array([0])).sum()) + 1
    padded = so3lr.jraph.pad_with_graphs(
        batch,
        n_node=int(batch.n_node.sum()) + 1,
        n_edge=int(batch.n_edge.sum()) + 1,
        n_graph=3,
        n_pairs=n_pairs,
    )
    inputs = so3lr.jraph_utils.graph_to_batch_fn(padded)
    for key in ('energy', 'forces'):
        inputs.pop(key, None)
    return inputs


def _warn_on_extreme_ratios(a0, c6):
    ratio_min, ratio_max = 0.05, 3.0
    combined = np.concatenate([c6, a0])
    if np.any((combined < ratio_min) | (combined > ratio_max)):
        warnings.warn(
            f'MBD-ML predicted alpha_0 or C6 ratios outside [{ratio_min}, '
            f'{ratio_max}]. Either the system is pathological or the model is '
            'extrapolating; the result may not be reliable.',
            stacklevel=3,
        )


def _ratios(atoms, model, jacobian=False):
    """Predict the ratios, and optionally their derivatives w.r.t. positions.

    The Jacobians are returned per Bohr, matching libMBD's gradients, while the
    model works in Angstrom.
    """
    so3lr = _so3lr()
    net, params, cutoff = _load_model(model)
    n_atoms = len(atoms)
    inputs = _graph_inputs(atoms, cutoff)
    positions = so3lr.jax.numpy.asarray(inputs['positions'])
    rest = {k: v for k, v in inputs.items() if k != 'positions'}

    def ratios_of(pos, other):
        out = net.apply(params, dict(positions=pos, **other))
        return out['hirshfeld_ratios'][:n_atoms], out['c6_ratios'][:n_atoms]

    a0, c6 = (np.asarray(x, dtype=float) for x in ratios_of(positions, rest))
    _warn_on_extreme_ratios(a0, c6)
    if not jacobian:
        return a0, c6, None, None
    da0, dc6 = (
        np.asarray(x, dtype=float)[:, :n_atoms, :] * Bohr
        for x in so3lr.jax.jacrev(ratios_of, argnums=0)(positions, rest)
    )
    return a0, c6, da0, dc6


def _shipped_models():
    """Names of the MBD-ML models so3lr ships, or () if so3lr is unavailable."""
    try:
        return _so3lr().MBD_ML_MODELS
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
        return _so3lr().model_path(model)
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
    a0, c6, _, _ = _ratios(atoms, model)
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


def mbd_properties_from_structure(
    atoms, beta, k_grid=None, model=DEFAULT_MODEL, ratio_response=True
):
    r"""Compute the MBD energy, forces and stress of a structure with MBD-ML ratios.

    :param atoms: ASE ``Atoms`` object
    :param float beta: MBD range-separation parameter
    :param k_grid: k-point grid, required for periodic systems
    :param model: name of an MBD-ML model shipped by so3lr, or the path of a
        directory holding one
    :param bool ratio_response: if True, add the term that accounts for the
        ratios themselves depending on the geometry (see below)

    Returns a dict with the energy ``'E'``, forces ``'F'`` and, for periodic
    systems, the stress ``'S'``, in atomic units.

    libMBD differentiates at fixed :math:`\alpha_0`, :math:`C_6` and
    :math:`R_\mathrm{vdw}`, but here those come from a model that depends on
    the geometry, so the bare gradient is not the gradient of the energy. The
    missing term,

    .. math::

        \sum_a
        \frac{\partial E}{\partial\alpha_a}\frac{\mathrm d\alpha_a}{\mathrm d\mathbf R}
        + \frac{\partial E}{\partial C_{6,a}}
          \frac{\mathrm dC_{6,a}}{\mathrm d\mathbf R}
        + \frac{\partial E}{\partial R_{\mathrm{vdw},a}}
          \frac{\mathrm dR_{\mathrm{vdw},a}}{\mathrm d\mathbf R},

    combines libMBD's vdW-parameter gradients with the model's own Jacobian,
    which jax supplies. It is not small: on water it is a tenth of the force,
    and without it a relaxation does not converge to a minimum of the energy it
    reports. ``ratio_response=False`` restores the older, inconsistent
    behaviour for comparison.

    The stress does not yet carry the corresponding lattice term, so for a
    periodic system ``'S'`` remains the fixed-ratio stress.
    """
    periodic = all(atoms.pbc)
    if any(atoms.pbc) and not periodic:
        raise ValueError(
            'MBD-ML supports fully periodic or fully non-periodic systems, '
            f'got pbc={tuple(bool(x) for x in atoms.pbc)}'
        )
    if periodic and k_grid is None:
        raise ValueError('k_grid must be given for periodic systems')

    a0_ratios, C6_ratios, da0_dR, dC6_dR = _ratios(
        atoms, model, jacobian=ratio_response
    )

    atom_pos = atoms.get_positions() / Bohr
    atom_species = atoms.get_chemical_symbols()
    n_atoms = len(atoms)
    lattice_vecs = atoms.cell[:, :] / Bohr if periodic else None

    a0_free, C6_free, _ = from_volumes(atom_species, np.ones(n_atoms))
    # the QDO vdW radius, which agrees with libMBD's R_vdw(TS) table to 1%
    Rvdw_free = 2.5 * a0_free ** (1 / 7)
    Rvdw = Rvdw_free * a0_ratios ** (1 / 3)

    a0 = a0_free * a0_ratios
    C6 = C6_free * C6_ratios

    confMBD = MBDGeom(coords=atom_pos, lattice=lattice_vecs, k_grid=k_grid)

    # mbd_energy() returns energies in Ha and energy gradients in Ha/Bohr
    results = confMBD.mbd_energy(
        a0,
        C6,
        R_vdw=Rvdw,
        beta=beta,
        damping='fermi,dip',
        variant='plain',
        force=True,
        vdw_params_grad=ratio_response,
    )
    if periodic:
        E, gradE, dE_dL = results[:3]
    else:
        (E, gradE), dE_dL = results[:2], None

    if ratio_response:
        dE_da0, dE_dC6, dE_dRvdw = results[-3:]
        # alpha_0 = alpha_free * v, C6 = C6_free * c, R_vdw = R_free * v**(1/3)
        dE_dv = dE_da0 * a0_free + dE_dRvdw * Rvdw_free * a0_ratios ** (-2 / 3) / 3
        gradE = (
            gradE
            + np.einsum('a,ajk->jk', dE_dv, da0_dR)
            + np.einsum('a,ajk->jk', dE_dC6 * C6_free, dC6_dR)
        )

    F = -gradE

    if periodic:
        S = compute_stress_from_lattice_gradient(lattice_vecs, atom_pos, dE_dL, -F)
        return {'E': E, 'F': F, 'S': S}
    return {'E': E, 'F': F}
