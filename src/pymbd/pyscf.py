# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
r"""Bridge between PySCF and pyMBD via Hirshfeld volume rescaling.

This optional module turns a converged PySCF Kohn--Sham object into an MBD (or
TS) dispersion energy. It computes per-atom Hirshfeld volume ratios

.. math::

    v_A = \frac{\int w_A(\mathbf r)\,|\mathbf r-\mathbf R_A|^3\,
                 \rho(\mathbf r)\,d\mathbf r}
               {\int |\mathbf r-\mathbf R_A|^3\,
                 \rho_A^\mathrm{free}(\mathbf r)\,d\mathbf r},
    \qquad
    w_A = \frac{\rho_A^\mathrm{free}}{\sum_B \rho_B^\mathrm{free}},

directly from the molecular density and a promolecule of spherically averaged
free-atom densities, and rescales the free-atom vdW parameters
(:math:`\alpha\propto v`, :math:`C_6\propto v^2`, :math:`R_\mathrm{vdw}\propto
v^{1/3}`) before calling :class:`pymbd.fortran.MBDGeom`.

The free-atom reference is a *spherically averaged Kohn--Sham* atom computed with
the same functional as the molecule (:class:`AtomSphAverageRKS`) in a large,
fixed basis (aug-cc-pVQZ by default). This matches the FHI-aims/libMBD
convention: for a one-electron atom (hydrogen) the reference is the diffuse,
self-consistent DFT density. Benchmarking against independent all-electron
FHI-aims MBD energies confirms this choice reproduces them best (MBD interaction
energies agree to ~0.02 kcal/mol). MBD is only meaningful on top of a
(semi)local/hybrid KS functional, so a Hartree--Fock mean-field is rejected.

Requires PySCF (``pip install pymbd[pyscf]``).
"""

from functools import lru_cache

import numpy as np
from pyscf import dft, gto
from pyscf.scf import atom_hf

from .fortran import MBDGeom
from .pymbd import vdw_params

__all__ = ['hirshfeld_volume_ratios', 'mbd_energy', 'AtomSphAverageRKS']

# MBD@rsSCS range-separation parameter beta per exchange-correlation functional.
_BETA_RSSCS = {'pbe': 0.83, 'pbe0': 0.85, 'hse': 0.85}

# Grid is processed in blocks: the promolecule AO array would otherwise be
# ngrid x nao of the (large) free-atom basis, which reaches several GB for a
# few dozen atoms in aug-cc-pVQZ.
_GRID_BLKSIZE = 8192
_AO_BLOCK_BYTES = 128 << 20


class AtomSphAverageRKS(dft.rks.KohnShamDFT, atom_hf.AtomSphAverageRHF):
    """Spherically averaged neutral free-atom solver with a Kohn--Sham potential.

    Composed the way PySCF builds ``RKS = KohnShamDFT + RHF``, but with the
    spherical-average ``eig``/``get_occ`` of :class:`pyscf.scf.atom_hf.
    AtomSphAverageRHF`, so the free-atom density is spherical (angular averaging
    of the Fock matrix plus fractional occupation of open shells) *and*
    consistent with the molecular functional.
    """

    def __init__(self, mol, xc='LDA,VWN'):
        atom_hf.AtomSphAverageRHF.__init__(self, mol)
        dft.rks.KohnShamDFT.__init__(self, xc)

    get_veff = dft.rks.get_veff
    energy_elec = dft.rks.energy_elec


@lru_cache(maxsize=None)
def _free_atom_dm(symb, basis, xc):
    """Spherically averaged free-atom density matrix for one element.

    Cached per ``(element, basis, xc)``, so each atomic SCF runs once per process.
    Since it only depends on the element and the basis, the result is expressed in
    exactly the AO ordering of that element's block of the promolecule (see
    :func:`_promolecule_mole`), where it is used. The one-electron atom (hydrogen)
    is kept as the self-consistent DFT density (diffuse, self-interaction and
    all), which is the MBD/TS convention. ``spin`` is set to the electron-count
    parity only to satisfy the Mole build; the spherical-average solver fixes the
    occupation by fractional filling.
    """
    mol = gto.M(
        atom=[[symb, (0.0, 0.0, 0.0)]],
        basis={symb: basis},
        spin=gto.charge(symb) % 2,
        verbose=0,
    )
    mf = AtomSphAverageRKS(mol, xc=xc)
    mf.kernel()
    return mf.make_rdm1()


def _promolecule_mole(mol, free_atom_basis):
    """Mole with the free-atom basis on every atom of ``mol`` (the promolecule).

    Only ever used to evaluate basis functions on a grid, never for an SCF, so
    ``spin`` just has to match the parity of the neutral-atom electron count.
    """
    return gto.M(
        atom=[
            [mol.atom_symbol(ia), tuple(mol.atom_coord(ia))] for ia in range(mol.natm)
        ],
        basis=free_atom_basis,
        unit='Bohr',
        spin=sum(gto.charge(mol.atom_symbol(ia)) for ia in range(mol.natm)) % 2,
        verbose=0,
    )


def hirshfeld_volume_ratios(
    mf, free_atom_xc='auto', free_atom_basis='aug-cc-pVQZ', grid_level=3
):
    """Per-atom Hirshfeld volume ratios from a converged PySCF Kohn--Sham object.

    :param mf: converged PySCF KS object (restricted or unrestricted)
    :param free_atom_xc: functional for the free-atom reference; ``'auto'`` (the
        default) uses the molecular functional, or pass an explicit functional
    :param str free_atom_basis: basis for the free-atom reference (fixed, large)
    :param int grid_level: fallback PySCF grid level, only used if ``mf`` has no
        integration grid of its own

    Returns an array of per-atom ratios of shape ``(natm,)``.

    Requires a KS-DFT mean-field; a Hartree--Fock ``mf`` (no exchange-correlation
    functional) is rejected. The molecular integration grid is taken from ``mf``.

    An isolated atom returns a ratio of exactly 1 only when ``free_atom_basis``
    matches the molecular basis (or both are at the basis-set limit); with the
    default fixed reference basis a small molecular basis leaves a systematic
    per-atom offset that largely cancels in interaction energies.
    """
    mol = mf.mol
    if free_atom_xc == 'auto':
        free_atom_xc = getattr(mf, 'xc', None)
    if not free_atom_xc:
        raise ValueError(
            'MBD Hirshfeld partitioning requires a KS-DFT mean-field; got one '
            'with no exchange-correlation functional (Hartree-Fock)'
        )
    ghosts = [ia for ia in range(mol.natm) if mol.atom_charge(ia) == 0]
    if ghosts:
        raise ValueError(
            'Hirshfeld partitioning is undefined for ghost atoms (no nucleus, '
            f'hence no free-atom reference); found at atom indices {ghosts}. '
            'Evaluate the dispersion on the ghost-free fragments: unlike the '
            'DFT interaction energy, the MBD term needs no counterpoise.'
        )
    grids = getattr(mf, 'grids', None)
    if grids is None or grids.coords is None:
        grids = dft.gen_grid.Grids(mol)
        grids.level = grid_level
        grids.build()
    dm = mf.make_rdm1()
    if np.ndim(dm) == 3:  # unrestricted -> total density
        dm = dm[0] + dm[1]

    # All free-atom densities come from one basis (the promolecule Mole) evaluated
    # in a single pass over the grid: the free-atom density matrices are
    # block-diagonal in it, so atom A's density is the contraction of A's own AO
    # block with its cached density matrix. Grid blocking keeps the AO array
    # bounded and lets the integrals be accumulated on the fly, so neither the
    # (natm, ngrid) promolecule matrix nor a full AO array is ever materialized.
    promol = _promolecule_mole(mol, free_atom_basis)
    aoslices = promol.aoslice_by_atom()
    free_dms = [
        _free_atom_dm(mol.atom_symbol(ia), free_atom_basis, free_atom_xc)
        for ia in range(mol.natm)
    ]
    atom_coords = mol.atom_coords()
    blksize = max(128, min(_GRID_BLKSIZE, _AO_BLOCK_BYTES // (8 * promol.nao)))

    v_eff = np.zeros(mol.natm)
    v_free = np.zeros(mol.natm)
    for i0 in range(0, len(grids.coords), blksize):
        coords = grids.coords[i0 : i0 + blksize]
        weights = grids.weights[i0 : i0 + blksize]
        ao_free = dft.numint.eval_ao(promol, coords)
        rho_free = np.empty((mol.natm, len(coords)))
        for ia in range(mol.natm):
            _, _, p0, p1 = aoslices[ia]
            ao_a = ao_free[:, p0:p1]
            rho_free[ia] = np.einsum('gi,ij,gj->g', ao_a, free_dms[ia], ao_a)
        np.clip(rho_free, 0, None, out=rho_free)
        rho_promol = rho_free.sum(axis=0)
        np.maximum(rho_promol, 1e-30, out=rho_promol)
        rho = dft.numint.eval_rho(mol, dft.numint.eval_ao(mol, coords), dm)
        for ia in range(mol.natm):
            r3 = np.linalg.norm(coords - atom_coords[ia], axis=1) ** 3
            w = rho_free[ia] / rho_promol  # Hirshfeld weight
            v_eff[ia] += np.einsum('g,g,g,g->', weights, w, r3, rho)
            v_free[ia] += np.einsum('g,g,g->', weights, r3, rho_free[ia])
    return v_eff / v_free


def mbd_energy(
    mf, beta=None, variant='rsscs', a=6.0, damping='fermi,dip', **ratio_kwargs
):
    """MBD dispersion energy (a.u.) for a converged PySCF Kohn--Sham object.

    :param mf: converged PySCF KS object
    :param float beta: MBD range-separation parameter; if ``None`` it is looked
        up from the molecular functional (currently PBE, PBE0, HSE)
    :param str variant: MBD variant passed to :meth:`~pymbd.fortran.MBDGeom.
        mbd_energy` (``'rsscs'`` or ``'plain'``)
    :param float a: MBD damping steepness
    :param str damping: MBD damping function
    :param ratio_kwargs: forwarded to :func:`hirshfeld_volume_ratios`

    Free-atom vdW parameters are taken from libMBD's built-in table and rescaled
    by the Hirshfeld volume ratios.
    """
    mol = mf.mol
    if beta is None:
        xc = getattr(mf, 'xc', '').lower().replace(' ', '')
        beta = _BETA_RSSCS.get(xc)
        if beta is None:
            raise ValueError(f'No default MBD beta for xc={xc!r}; pass beta explicitly')
    ratios = hirshfeld_volume_ratios(mf, **ratio_kwargs)
    species = [mol.atom_symbol(i) for i in range(mol.natm)]
    alpha_0 = np.array([vdw_params[s]['alpha_0(TS)'] for s in species]) * ratios
    C6 = np.array([vdw_params[s]['C6(TS)'] for s in species]) * ratios**2
    R_vdw = np.array([vdw_params[s]['R_vdw(TS)'] for s in species]) * ratios ** (1 / 3)
    return MBDGeom(mol.atom_coords()).mbd_energy(
        alpha_0, C6, R_vdw, beta=beta, a=a, variant=variant, damping=damping
    )
