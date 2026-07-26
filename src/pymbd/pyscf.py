# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
r"""Bridge between PySCF and pyMBD via Hirshfeld volume rescaling.

This optional module turns a converged PySCF mean-field object into an MBD (or
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
self-consistent DFT density, *not* the exact density -- benchmarking against
independent all-electron FHI-aims MBD energies confirms this choice reproduces
them best (MBD interaction energies agree to ~0.02 kcal/mol).

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
        self.init_guess = 'vsap'

    get_veff = dft.rks.get_veff
    energy_elec = dft.rks.energy_elec


@lru_cache(maxsize=None)
def _free_atom_dm(symb, basis, xc, h_exact=False):
    """Spherically averaged free-atom density matrix for one element (cached).

    ``xc=None`` gives a Hartree--Fock reference; otherwise a Kohn--Sham reference
    with functional ``xc``. A one-electron atom is routed to the exact one-body
    solver only in the HF path or when ``h_exact`` is set; in the KS path it is
    kept as the (diffuse, self-interaction-carrying) self-consistent DFT density,
    which is the MBD/TS convention.
    """
    atm = gto.M(
        atom=[[symb, (0.0, 0.0, 0.0)]],
        basis={symb: basis},
        spin=gto.charge(symb) % 2,
        verbose=0,
    )
    if atm.nelectron == 1 and (xc is None or h_exact):
        mf = atom_hf.AtomHF1e(atm)
    elif xc is None:
        mf = atom_hf.AtomSphAverageRHF(atm)
    else:
        mf = AtomSphAverageRKS(atm, xc=xc)
    mf.kernel()
    return mf.make_rdm1()


def _free_atom_density_on_grid(mol, coords, xc, free_atom_basis, h_exact):
    """Return per-atom free-atom densities on ``coords`` and their sum."""
    rho_free = np.zeros((mol.natm, len(coords)))
    for ia in range(mol.natm):
        symb = mol.atom_symbol(ia)
        dm = _free_atom_dm(symb, free_atom_basis, xc, h_exact)
        atm = gto.M(
            atom=[[symb, mol.atom_coord(ia, unit='Bohr')]],
            basis={symb: free_atom_basis},
            unit='Bohr',
            spin=gto.charge(symb) % 2,
            verbose=0,
        )
        ao = dft.numint.eval_ao(atm, coords)
        rho_free[ia] = dft.numint.eval_rho(atm, ao, dm)
    rho_free = np.clip(rho_free, 0, None)
    return rho_free, rho_free.sum(axis=0)


def hirshfeld_volume_ratios(
    mf,
    grid_level=3,
    free_atom_xc='auto',
    free_atom_basis='aug-cc-pVQZ',
    h_exact=False,
):
    """Per-atom Hirshfeld volume ratios from a converged PySCF mean-field.

    :param mf: converged PySCF SCF/KS object (restricted or unrestricted)
    :param int grid_level: PySCF DFT grid level for the integration
    :param free_atom_xc: functional for the free-atom reference; ``'auto'`` uses
        the molecular functional (HF if ``mf`` has none), ``None`` forces HF, or
        an explicit functional string
    :param str free_atom_basis: basis for the free-atom reference (fixed, large)
    :param bool h_exact: use the exact one-electron density for hydrogen instead
        of the self-consistent DFT one (off by default; see module docstring)

    Returns an array of per-atom ratios of shape ``(natm,)``.

    An isolated atom returns a ratio of exactly 1 only when ``free_atom_basis``
    matches the molecular basis (or both are at the basis-set limit); with the
    default fixed reference basis a small molecular basis leaves a systematic
    per-atom offset that largely cancels in interaction energies.
    """
    mol = mf.mol
    if free_atom_xc == 'auto':
        free_atom_xc = getattr(mf, 'xc', None)
    grids = dft.gen_grid.Grids(mol)
    grids.level = grid_level
    grids.build()
    coords, weights = grids.coords, grids.weights

    dm = mf.make_rdm1()
    if np.ndim(dm) == 3:  # unrestricted -> total density
        dm = dm[0] + dm[1]
    ao = dft.numint.eval_ao(mol, coords)
    rho = dft.numint.eval_rho(mol, ao, dm)

    rho_free, rho_promol = _free_atom_density_on_grid(
        mol, coords, free_atom_xc, free_atom_basis, h_exact
    )
    rho_promol = np.where(rho_promol > 1e-30, rho_promol, 1e-30)

    atom_coords = mol.atom_coords()
    ratios = np.empty(mol.natm)
    for ia in range(mol.natm):
        r3 = np.linalg.norm(coords - atom_coords[ia], axis=1) ** 3
        w = rho_free[ia] / rho_promol
        v_eff = np.einsum('g,g,g,g->', weights, w, r3, rho)
        v_free = np.einsum('g,g,g->', weights, r3, rho_free[ia])
        ratios[ia] = v_eff / v_free
    return ratios


def mbd_energy(
    mf, beta=None, variant='rsscs', a=6.0, damping='fermi,dip', **ratio_kwargs
):
    """MBD dispersion energy (a.u.) for a converged PySCF mean-field.

    :param mf: converged PySCF SCF/KS object
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
        xc = getattr(mf, 'xc', 'hf').lower().replace(' ', '')
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
