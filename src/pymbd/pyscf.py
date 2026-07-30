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

The free-atom reference is a *spherically averaged Kohn--Sham* atom
(:class:`AtomSphAverageRKS`) computed with the same functional *and the same basis*
as the molecule. Sharing the basis is what makes the ratio a pure measure of the
chemical environment: the reference then cancels the basis-set incompleteness that
is also present in the molecular density, so well-separated atoms come out at
exactly 1 in any basis. For a one-electron atom (hydrogen) the reference is the
diffuse, self-consistent DFT density, matching the FHI-aims/libMBD convention.
Benchmarking against independent all-electron FHI-aims MBD energies confirms these
choices: MBD interaction energies agree to ~0.003 kcal/mol (against ~0.016 with a
fixed aug-cc-pVQZ reference). MBD is only meaningful on top of a
(semi)local/hybrid KS functional, so a Hartree--Fock mean-field is rejected.

Requires PySCF (``pip install pymbd[pyscf]``)::

    from pyscf import dft, gto

    from pymbd.pyscf import mbd_energy

    mol = gto.M(atom='He 0 0 0; He 0 0 3', basis='def2-tzvp')
    mf = dft.RKS(mol, xc='PBE').run()
    ene = mbd_energy(mf)  # beta is taken from the functional
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

# Grid is processed in blocks: a whole-grid AO array is ngrid x nao, which runs
# into several GB for a few dozen atoms in a triple-zeta or larger basis.
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


def _hashable_basis(basis):
    """Make a basis spec hashable (nested lists -> tuples) for use as a cache key.

    A parsed per-element basis (``mol._basis[symb]``) is a nested list; PySCF
    accepts the tuple form just as well, so the result doubles as the spec passed
    back to :func:`pyscf.gto.M`.
    """
    if isinstance(basis, (list, tuple)):
        return tuple(_hashable_basis(b) for b in basis)
    return basis  # str basis name, or a number inside a parsed basis


@lru_cache(maxsize=None)
def _free_atom_dm(symb, basis, xc):
    """Spherically averaged free-atom density matrix for one element.

    Cached per ``(element, basis, xc)``, so each atomic SCF runs once per process.
    Because it depends only on the element and the basis, the result comes out in
    exactly the AO ordering of that element's block of the molecule it is used
    with, which is what lets it be contracted against a slice of the molecular AOs.
    The one-electron atom (hydrogen) is kept as the self-consistent DFT density
    (diffuse, self-interaction and all), which is the MBD/TS convention. ``spin``
    is set to the electron-count parity only to satisfy the Mole build; the
    spherical-average solver fixes the occupation by fractional filling.
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


def hirshfeld_volume_ratios(mf):
    """Per-atom Hirshfeld volume ratios from a converged PySCF Kohn--Sham object.

    :param mf: converged PySCF KS object (restricted or unrestricted)

    Returns an array of per-atom ratios of shape ``(natm,)``.

    Requires an all-electron KS-DFT mean-field; Hartree--Fock (no
    exchange-correlation functional), ghost atoms and ECPs are rejected. The
    functional, basis and integration grid are all taken from ``mf``.

    The free-atom reference is deliberately computed in the *molecular* basis. Any
    other choice mixes basis-set incompleteness into the ratio: a well-separated
    argon dimer in def2-SVP gives 0.89 against a fixed aug-cc-pVQZ reference (a
    21% error in :math:`C_6`) but 1.000000 against a def2-SVP one, since the
    reference then cancels the incompleteness that is also in the molecular
    density. The two agree only at the basis-set limit.
    """
    mol = mf.mol
    xc = getattr(mf, 'xc', None)
    if not xc:
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
    if mol.has_ecp():
        raise ValueError(
            'Hirshfeld volumes need the all-electron density, but this molecule '
            'uses an ECP, so its core density is missing and cannot be matched '
            'by a neutral free-atom reference'
        )
    grids = mf.grids  # a KS object always carries one
    if grids.coords is None:
        grids.build()
    dm = mf.make_rdm1()
    if np.ndim(dm) == 3:  # unrestricted -> total density
        dm = dm[0] + dm[1]

    # Since the free-atom reference uses the molecular basis, the promolecule is
    # just mol and one AO evaluation serves both the free-atom densities and the
    # molecular one: the free-atom density matrices are block-diagonal in that
    # basis, so atom A's density is the contraction of A's own AO block with its
    # cached density matrix. Grid blocking keeps the AO array bounded and lets the
    # integrals be accumulated on the fly, so neither a whole-grid AO array nor an
    # (natm, ngrid) promolecule matrix is ever materialized.
    symbols = [mol.atom_symbol(ia) for ia in range(mol.natm)]
    free_dms = [
        _free_atom_dm(symb, _hashable_basis(mol._basis[symb]), xc) for symb in symbols
    ]
    aoslices = mol.aoslice_by_atom()
    atom_coords = mol.atom_coords()
    blksize = max(128, min(_GRID_BLKSIZE, _AO_BLOCK_BYTES // (8 * mol.nao)))

    v_eff = np.zeros(mol.natm)
    v_free = np.zeros(mol.natm)
    for i0 in range(0, len(grids.coords), blksize):
        coords = grids.coords[i0 : i0 + blksize]
        weights = grids.weights[i0 : i0 + blksize]
        ao = dft.numint.eval_ao(mol, coords)
        rho_free = np.empty((mol.natm, len(coords)))
        for ia in range(mol.natm):
            _, _, p0, p1 = aoslices[ia]
            ao_a = ao[:, p0:p1]
            # contract via a matrix product, as eval_rho does internally: going
            # through BLAS is ~10x faster than a three-operand einsum here, and
            # the block structure avoids the (much larger) dense nao x nao product
            rho_free[ia] = np.einsum('gi,gi->g', ao_a.dot(free_dms[ia]), ao_a)
        np.clip(rho_free, 0, None, out=rho_free)
        rho_promol = rho_free.sum(axis=0)
        np.maximum(rho_promol, 1e-30, out=rho_promol)
        rho = dft.numint.eval_rho(mol, ao, dm)
        for ia in range(mol.natm):
            r3 = np.linalg.norm(coords - atom_coords[ia], axis=1) ** 3
            w = rho_free[ia] / rho_promol  # Hirshfeld weight
            v_eff[ia] += np.einsum('g,g,g,g->', weights, w, r3, rho)
            v_free[ia] += np.einsum('g,g,g->', weights, r3, rho_free[ia])
    return v_eff / v_free


def mbd_energy(mf, beta=None, variant='rsscs', a=6.0, damping='fermi,dip'):
    """MBD dispersion energy (a.u.) for a converged PySCF Kohn--Sham object.

    :param mf: converged PySCF KS object
    :param float beta: MBD range-separation parameter; if ``None`` it is looked
        up from the molecular functional (currently PBE, PBE0, HSE)
    :param str variant: MBD variant passed to :meth:`~pymbd.fortran.MBDGeom.
        mbd_energy` (``'rsscs'`` or ``'plain'``)
    :param float a: MBD damping steepness
    :param str damping: MBD damping function

    Free-atom vdW parameters are taken from libMBD's built-in table and rescaled
    by the Hirshfeld volume ratios.
    """
    mol = mf.mol
    if beta is None:
        xc = getattr(mf, 'xc', '').lower().replace(' ', '')
        beta = _BETA_RSSCS.get(xc)
        if beta is None:
            raise ValueError(f'No default MBD beta for xc={xc!r}; pass beta explicitly')
    ratios = hirshfeld_volume_ratios(mf)
    species = [mol.atom_symbol(i) for i in range(mol.natm)]
    alpha_0 = np.array([vdw_params[s]['alpha_0(TS)'] for s in species]) * ratios
    C6 = np.array([vdw_params[s]['C6(TS)'] for s in species]) * ratios**2
    R_vdw = np.array([vdw_params[s]['R_vdw(TS)'] for s in species]) * ratios ** (1 / 3)
    return MBDGeom(mol.atom_coords()).mbd_energy(
        alpha_0, C6, R_vdw, beta=beta, a=a, variant=variant, damping=damping
    )
