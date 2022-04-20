import os

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


class DispersionCorrectionCalculator(Calculator):
    """ASE's Calculator to combine two calculators in a simple additive scheme.

    This :class:`~ase.calculators.calculator.Calculator` hosts several
    calculators and simply adds energy, force, stress contributions. Has special
    features for intercalc communication for dispersion corrections.

    TODO: short-range has to implement hirsh_volrat property
    """

    implemented_properties = ['energy', 'forces', 'stress']

    valid_args = (
        'qm_calculator',  # array of calculators
        'mm_calculator',  # array of calculators
        'hirbulk',  # float
        'hirlast',
    )  # array of tuples

    def __init__(
        self, restart=None, label=os.curdir, atoms=None, reset=False, **kwargs
    ):
        """Assign individual calculators.

        Parameters
        ==========
        calcs:        list
            list of attached calculators

        qm_calculator:     list of members of a Class defining a Calculator
            ase qm-calculator for each qm region

        mm_calculator:      member of a Class defining a Calculator
            ase mm-calculator for the mm region (the whole system)

        hirbulk:        float
            Default value for Hirshfeld volume ratio if mm_calc can take
            hirshfeld ratios for vdW calculations
        """

        # Set any keyword arguments
        for arg, val in kwargs.items():
            if arg in self.valid_args:
                setattr(self, arg, val)
            else:
                raise RuntimeError(
                    'unknown keyword arg "%s" : not in %s' % (arg, self.valid_args)
                )

        # Check the user input for missing information
        error_head = ' +----------------** QMME WARNING **----------------+'
        error_tail = ' +--------------------------------------------------+'
        for arg in self.valid_args:
            if arg not in kwargs.keys():
                if arg == 'qm_calculator':
                    print(error_head)
                    print(' |  Keyword  qm_calculator  not specified.')
                    print(error_tail + '\n')
                if arg == 'mm_calculator':
                    print(error_head)
                    print(' |  Keyword  mm_calculator  not specified.')
                    print(error_tail + '\n')

        # reset to v_hirsh/v_free = 1 if not give
        if not hasattr(self, 'hirbulk'):
            self.hirbulk = 1.0

        # Check if the QM-calculator support Hirshfeld-paritioning and inform
        # user about missing Hirshfeld-Treatment.
        if str(self.mm_calculator) == 'MBD':
            qmcalc = self.qm_calculator
            if not hasattr(qmcalc, 'get_hirsh_volrat'):
                print(error_head)
                print(
                    ' |  QM calculator '
                    + str(qmcalc.__class__.__name__)
                    + ' does not support Hirshfeld-'
                )
                print(' |  partitioning. Defaulting to v_hirsh/v_free = bulk_value')
                print(error_tail + '\n')

        if str(self.qm_calculator) == 'Aims':
            self.reset = True
        if str(self.qm_calculator) in ['SpkVdwCalculator', 'SpkCalculator']:
            self.reset = False

        # Initialization of empty data structures of calculation results.
        # This is required for some internal reconstruction purposes and
        # error processing.
        self.nopt = 0

        Calculator.__init__(self, restart, label, **kwargs)

    def calculate(self, atoms, properties=('energy',), system_changes=all_changes):
        """Do the calculation."""
        self.atoms = atoms.copy()

        if self.calculation_required(atoms, properties):
            Calculator.calculate(self, atoms)
            qm_atoms = self.atoms.copy()
            mm_atoms = self.atoms.copy()

            del qm_atoms.calc
            del mm_atoms.calc
            # QM part
            if self.reset:
                self.qm_calculator.reset()
            qm_atoms.set_calculator(self.qm_calculator)

            qm_energy = qm_atoms.get_potential_energy()
            if 'forces' in properties:
                qm_forces = qm_atoms.get_forces()
            if 'stress' in properties and all(qm_atoms.get_pbc()):
                qm_stress = qm_atoms.get_stress()

            if hasattr(self.qm_calculator, 'get_hirsh_volrat'):
                self.hirlast = qm_atoms.calc.get_hirsh_volrat()
            else:
                self.hirlast = np.ones(len(atoms)) * self.hirbulk

            # MM part
            mm_atoms.set_calculator(self.mm_calculator)
            if hasattr(self.mm_calculator, 'set_hirshfeld'):
                mm_atoms.calc.set_hirshfeld(self.hirlast)

            mm_energy = mm_atoms.get_potential_energy()
            if 'forces' in properties:
                mm_forces = mm_atoms.get_forces()
            if 'stress' in properties and all(qm_atoms.get_pbc()):
                mm_stress = mm_atoms.get_stress()

            self.results['energy'] = self.results['free_energy'] = qm_energy + mm_energy
            if 'forces' in properties:
                self.results['forces'] = qm_forces + mm_forces
            if 'stress' in properties and all(self.atoms.get_pbc()):
                self.results['stress'] = qm_stress + mm_stress

            self.nopt += 1
