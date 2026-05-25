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
        n_freq=15,
        beta=0.83,  # PBE default value
        ts_sr=0.94,  # PBE default value
        do_rpa=False,
        **kwargs,
    ):
        scheme = scheme.upper()
        params = params.upper()
        assert scheme in ['VDW', 'MBD']
        self.k_grid = k_grid
        self.scheme = scheme
        self.params = params
        self.n_freq = n_freq
        self.beta = beta
        self.ts_sr = ts_sr
        self.do_rpa = do_rpa
        if self.params == 'TSSURF':
            self.params = 'TSsurf'
        self.hirsh_volrat = None
        super().__init__(restart, label, atoms, **kwargs)

    def calculate(self, atoms, properties=('energy',), system_changes=all_changes):
        """Do the calculation."""
        if not self.calculation_required(atoms, properties):
            return
        super().calculate(atoms, properties, system_changes)
        do_force = 'forces' in properties or 'stress' in properties
        hirsh_volrat = (
            self.hirsh_volrat
            if self.hirsh_volrat is not None
            else np.ones(len(atoms), dtype=np.float)
        )
        alpha_0, C6, R_vdw = from_volumes(
            atoms.get_chemical_symbols(), hirsh_volrat, kind=self.params
        )
        geom = MBDGeom(
            coords=atoms.positions / units.Bohr,
            lattice=(
                np.array(atoms.get_cell()) / units.Bohr
                if all(atoms.get_pbc())
                else None
            ),
            k_grid=self.k_grid,
            n_freq=self.n_freq,
            do_rpa=self.do_rpa,
        )
        if self.scheme == 'MBD':
            result = geom.mbd_energy(alpha_0, C6, R_vdw, beta=self.beta, force=do_force)
        elif self.scheme == 'VDW':
            result = geom.ts_energy(
                alpha_0, C6, R_vdw, sR=self.ts_sr, d=20.0, force=do_force
            )
        self.results['energy'] = self.results['free_energy'] = (
            result if not do_force else result[0]
        ) * units.Hartree
        if do_force:
            gradients = result[1] * units.Hartree / units.Bohr
            self.results['forces'] = -gradients
            if 'stress' in properties and all(self.atoms.get_pbc()):
                latt_gradients = result[2] * units.Hartree / units.Bohr
                stress = np.dot(atoms.get_cell(), latt_gradients.transpose()) + np.dot(
                    atoms.get_positions().transpose(), gradients
                )
                stress = stress / (atoms.get_volume())
                self.results['stress'] = stress

    def set_hirshfeld(self, hirsh_volrat):
        """Set Hirshfeld volumes."""
        self.hirsh_volrat = hirsh_volrat
