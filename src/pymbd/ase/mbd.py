import os

import ase.units as units
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from ..fortran import MBDGeom
from ..pymbd import from_volumes


class MBD(Calculator):
    """ASE's :class:`~ase.calculators.calculator.Calculator` for Libmbd."""

    name = 'MBD'
    implemented_properties = ['energy', 'forces', 'stress', 'free_energy']

    def __init__(
        self,
        restart=None,
        label=os.curdir,
        atoms=None,
        k_grid=None,
        scheme='VDW',  # 'VDWw' or 'MBD', default is VDW'
        params='TS',  # either TS or TSSURF
        nfreq=15,
        beta=0.83,  # PBE default value
        ts_sr=0.94,  # PBE default value
        do_rpa=False,
        **kwargs,
    ):
        self.k_grid = k_grid
        self.scheme = scheme.upper()
        self.params = params.upper()
        self.nfreq = nfreq
        self.beta = beta
        self.ts_sr = ts_sr
        self.do_rpa = do_rpa
        if self.params == 'TSSURF':
            self.params = 'TSsurf'
        self.hirshvolrat_is_set = False
        self.hirsh_volrat = None
        Calculator.__init__(self, restart, label, atoms, **kwargs)

    def calculate(self, atoms, properties=('energy',), system_changes=all_changes):
        """Do the calculation."""
        self.atoms = atoms.copy()

        if all(atoms.get_pbc()):
            lattice = np.array(atoms.get_cell()) / units.Bohr
        else:
            lattice = None

        if not self.hirshvolrat_is_set:
            self.hirsh_volrat = np.ones(len(atoms), dtype=np.float)

        if self.calculation_required(atoms, properties):
            Calculator.calculate(self, atoms)

            self.alpha_0, self.C6, self.R_vdw = from_volumes(
                atoms.get_chemical_symbols(), self.hirsh_volrat, kind=self.params
            )

            if 'forces' in properties or 'stress' in properties:
                do_force = True
            else:
                do_force = False

            mbdgeom = MBDGeom(
                coords=atoms.positions / units.Bohr,
                lattice=lattice,
                k_grid=self.k_grid,
                n_freq=self.nfreq,
                do_rpa=self.do_rpa,
            )

            if self.scheme == 'MBD':
                energy = mbdgeom.mbd_energy(
                    self.alpha_0, self.C6, self.R_vdw, beta=self.beta, force=do_force
                )
            elif self.scheme == 'VDW':
                energy = mbdgeom.ts_energy(
                    self.alpha_0,
                    self.C6,
                    self.R_vdw,
                    sR=self.ts_sr,
                    d=20.0,
                    force=do_force,
                )
            else:
                raise ValueError('mbd: scheme needs to be MBD or VDW')

            if not do_force:
                self.results['energy'] = energy * units.Hartree
            else:
                self.results['energy'] = energy[0] * units.Hartree
                gradients = energy[1] * units.Hartree / units.Bohr
                self.results['forces'] = -gradients
                if 'stress' in properties and all(self.atoms.get_pbc()):
                    lattgradients = energy[2] * units.Hartree / units.Bohr
                    stress = np.dot(
                        atoms.get_cell(),
                        lattgradients.transpose(),
                    ) + np.dot(atoms.get_positions().transpose(), gradients)
                    stress = stress / (atoms.get_volume())
                    self.results['stress'] = stress
            self.results['free_energy'] = self.results['energy']

    def set_hirshfeld(self, hirsh_volrat):
        """Set Hirshfeld volumes."""
        self.hirshvolrat_is_set = True
        self.hirsh_volrat = hirsh_volrat
