.. _examples:

=======================
Pymbd with ASE and GPAW
=======================

Example of how to calculate the polarizability of a nitrogen molecule
and how to calculate many-body dispersion corrections with Pymbd.


.. code:: python
from ase.build import molecule
from ase.units import Bohr, Ha
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
