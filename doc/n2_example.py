from ase.build import molecule
from ase.units import Bohr, Ha
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac
from gpaw.cluster import Cluster
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from pymbd import molecular_polarizability
from pymbd.fortran import MBDGeom

# make ASE atoms object
atoms = Cluster(molecule('N2'))

# initialize GPAW calculator
h = 0.2  # set grid spacing for GPAW calculation
atoms.minimal_box(4., h=h)
c = GPAW(xc='PBE', h=h, nbands=-6, occupations=FermiDirac(width=0.1))
atoms.calc = c

# get the ground state energy and forces
energy = atoms.get_potential_energy()
forces = atoms.get_forces()

# calculate the Hirshfeld volumes
hf = HirshfeldPartitioning(atoms.calc)
volumes = hf.get_effective_volume_ratios()

# convert input paramenters from ASE to pymbd
coords = atoms.positions / Bohr
species = atoms.get_chemical_symbols()
beta = 0.83  # default for PBE xc-functional

# static polarizability tensor for nitrogen in Bohr**3
alpha = molecular_polarizability(coords, species, volumes, beta)

# dispersion correction for energy and forces
mbd_atoms = MBDGeom(coords, lattice=None)
mbd_energy, mbd_forces = mbd_atoms.mbd_energy_species(species,
                                                      volumes,
                                                      beta,
                                                      force=True)

# convert from Hartree to eV
mbd_energy *= Ha
# convert from Ha / Bohr to eV/ Angstrom
mbd_forces *= Ha / Bohr

# add dispersion correction to ground state energy to get total energy and forces
total_energy = energy + mbd_energy
total_forces = forces + mbd_forces

f = paropen('n2_results', 'a')

print('Results for nitrogen molecule', file=f)
print('Molecular polarizability [Bohr**3]', file=f)
print(alpha, file=f)

print('Energy without dispersion correction [eV]', file=f)
print(energy, file=f)
print('Dispersion energy [eV]', file=f)
print(mbd_energy, file=f)
print('Total energy with dispersion correction [eV]', file=f)
print(total_energy, file=f)

print('Forces without dispersion correction [eV/Ang]', file=f)
print(forces, file=f)
print('Dispersion correction [eV/Ang]', file=f)
print(mbd_forces, file=f)
print('Total forces with dispersion correction [eV/Ang]', file=f)
print(total_forces, file=f)

f.close()

