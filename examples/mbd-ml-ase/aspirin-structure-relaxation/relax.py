#!/usr/bin/env python3

import ase.calculators.mixing
import jax
import numpy as np
from ase.calculators.aims import Aims
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write
from ase.optimize import BFGS

# For pymbd
from ase.units import Bohr, Hartree

import pymbd.mbd_ml

jax.config.update('jax_platform_name', 'cpu')


def get_deltak_for_geo(atoms, delta_k):
    recip_latt_vecs = atoms.get_reciprocal_cell()
    recip_latt_norms = np.linalg.norm(recip_latt_vecs, axis=1)
    k_grid = np.rint(recip_latt_norms / delta_k).astype(np.int64)

    return tuple([int(k) for k in k_grid])


class C6So3lrCalculator(Calculator):
    """So3lrCalculator that outputs MBD-correction of energy and forces"""

    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, k_grid=None, **kwargs):
        super().__init__(**kwargs)
        self.k_grid = k_grid
        self.is_mbd_computed = False
        # Total cutoff of c6-so3lr model in Angstrom (hardcoded for now)
        self.total_cutoff_hardcoded = 4.0

    def _compute_mbd_correction_from_ratios(
        self, atoms, k_grid, force_consistent=False
    ):

        prop_dict = pymbd.mbd_ml.mbd_properties_from_structure(
            atoms, beta=0.81, k_grid=k_grid
        )
        # libmbd outputs in Hartree and Bohr units -> convert to eV and A
        E = prop_dict['E'] * Hartree
        F = prop_dict['F'] * Hartree / Bohr
        S = prop_dict['S'] * Hartree / (Bohr**3)

        self.is_mbd_computed = True

        self.results['energy'] = float(E)
        self.results['forces'] = np.array(F)
        self.results['stress'] = np.array(S)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        properties = properties or ['energy', 'forces']
        super().calculate(atoms, properties, system_changes)
        self._compute_mbd_correction_from_ratios(atoms, self.k_grid)

    # Override get_potential_energy function
    def get_potential_energy(self, atoms=None, force_consistent=False):
        if atoms is None:
            atoms = self.atoms

        if not self.is_mbd_computed:
            self._compute_mbd_correction_from_ratios(
                atoms, self.k_grid, force_consistent
            )

        # get_potential_energy() now returns the MBD energy correction
        return self.results['energy']

    # Override get_forces
    def get_forces(self, atoms=None):
        if atoms is None:
            atoms = self.atoms

        if not self.is_mbd_computed:
            self._compute_mbd_correction_from_ratios(atoms, self.k_grid)

        # get_forces() now returns the MBD atomic force correction
        return self.results['forces']

    def get_stress(self, atoms=None):
        if atoms is None:
            atoms = self.atoms

        if not self.is_mbd_computed:
            self._compute_mbd_correction_from_ratios(atoms, self.k_grid)

        # get_stress() now returns the MBD stress tensor correction (3x3)
        return self.results['stress']

    def check_state(self, atoms, tol=1e-15):
        system_changes = super().check_state(atoms, tol)

        if system_changes:
            # Make sure that mbd correction is computed whenenever the structure changes
            self.is_mbd_computed = False

        return system_changes


atoms = read('geometry_original.in', format='aims')  # Replace with your structure file
k_grid = get_deltak_for_geo(atoms, 0.03)

# Create combined calculator between Aims() and C6So3lrCalculator(), where Aims()
# computes the total DFT energy and DFT atomic forces and C6So3lrCalculator()
# computes the MBD correction for energy and forces.

# Set up FHI-aims calculator
calc_aims = Aims(
    xc='pbe',
    spin='none',
    relativistic=('atomic_zora', 'scalar'),
    occupation_type=('gaussian', 0.05),
    basis_threshold=1e-5,
    k_grid=k_grid,
    species_dir='/path/to/fhi-aims/basisset/tight',
)

calc_c6so3lr = C6So3lrCalculator(k_grid=k_grid)

calc_mixed = ase.calculators.mixing.SumCalculator([calc_aims, calc_c6so3lr])

# Attach DFT+So3lr calculator to atoms
atoms.calc = calc_mixed

optimizer = BFGS(
    atoms,
    trajectory='relaxation.traj',  # Save trajectory
    logfile='optimization.log',
)

# Perform geometry optimization
print('Starting geometry relaxation...')
print(f'Initial energy: {atoms.get_potential_energy():.6f} eV')

# Run optimization (adjust fmax for desired convergence)
optimizer.run(fmax=0.005)  # Force convergence criterion in eV/Å

print('Geometry relaxation completed!')
print(f'Final energy: {atoms.get_potential_energy():.6f} eV')

# Save optimized structure
write('optimized_structure.extxyz', atoms)
write('optimized_structure.in', atoms, format='aims')

print(
    "Optimized structure saved as 'optimized_structure.extxyz' and 'optimized_structure.in'"
)
