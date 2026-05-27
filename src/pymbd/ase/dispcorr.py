import os
import warnings

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


class DispersionCorrectionCalculator(Calculator):
    """ASE's Calculator to combine two calculators in a simple additive scheme.

    This :class:`~ase.calculators.calculator.Calculator` hosts two calculators
    (a "QM" base calculator and a dispersion-correction "MM" calculator) and
    simply adds their energy, force and stress contributions. In addition, if
    the QM calculator exposes Hirshfeld volume ratios through a
    ``get_hirsh_volrat`` method, and the MM calculator accepts them through a
    ``set_hirshfeld`` method (as :class:`~pymbd.ase.MBD` does), the ratios are
    forwarded from the QM to the MM calculator before each evaluation. If the
    QM calculator does not provide Hirshfeld ratios, the value of ``hirbulk``
    (a float or an array of length ``len(atoms)``) is used as a fallback.
    """

    implemented_properties = ['energy', 'free_energy', 'forces', 'stress']

    def __init__(
        self,
        qm_calculator,
        mm_calculator,
        hirbulk=1.0,
        restart=None,
        label=os.curdir,
        atoms=None,
        **kwargs,
    ):
        """Initialize the calculator.

        :param qm_calculator: ASE calculator providing the base energy/forces.
            May expose Hirshfeld volume ratios via a ``get_hirsh_volrat``
            method, which will be forwarded to ``mm_calculator``.
        :param mm_calculator: ASE calculator providing the dispersion
            correction. May expose a ``set_hirshfeld`` method to receive the
            Hirshfeld volume ratios.
        :param hirbulk: Default Hirshfeld volume ratio (or array of ratios of
            length ``len(atoms)``) used when ``qm_calculator`` does not provide
            ``get_hirsh_volrat``. Defaults to 1.0.
        """
        self.qm_calculator = qm_calculator
        self.mm_calculator = mm_calculator
        self.hirbulk = hirbulk
        self.hirlast = None
        if getattr(self.mm_calculator, 'name', None) == 'MBD' and not hasattr(
            self.qm_calculator, 'get_hirsh_volrat'
        ):
            warnings.warn(
                'QM calculator {} does not provide Hirshfeld volume ratios; '
                'falling back to hirbulk={}.'.format(
                    self.qm_calculator.__class__.__name__, self.hirbulk
                ),
                stacklevel=2,
            )
        Calculator.__init__(self, restart=restart, label=label, atoms=atoms, **kwargs)

    def calculate(self, atoms, properties=('energy',), system_changes=all_changes):
        """Do the calculation."""
        Calculator.calculate(self, atoms, properties, system_changes)
        qm_atoms = self.atoms.copy()
        mm_atoms = self.atoms.copy()
        want_forces = 'forces' in properties
        want_stress = 'stress' in properties and all(self.atoms.get_pbc())

        # QM part
        qm_atoms.calc = self.qm_calculator
        qm_energy = qm_atoms.get_potential_energy()
        qm_forces = qm_atoms.get_forces() if want_forces else None
        qm_stress = qm_atoms.get_stress() if want_stress else None

        if hasattr(self.qm_calculator, 'get_hirsh_volrat'):
            self.hirlast = self.qm_calculator.get_hirsh_volrat()
        else:
            self.hirlast = np.ones(len(atoms)) * self.hirbulk

        # MM part
        mm_atoms.calc = self.mm_calculator
        if hasattr(self.mm_calculator, 'set_hirshfeld'):
            self.mm_calculator.set_hirshfeld(self.hirlast)

        mm_energy = mm_atoms.get_potential_energy()
        mm_forces = mm_atoms.get_forces() if want_forces else None
        mm_stress = mm_atoms.get_stress() if want_stress else None

        self.results['energy'] = self.results['free_energy'] = qm_energy + mm_energy
        if want_forces:
            self.results['forces'] = qm_forces + mm_forces
        if want_stress:
            self.results['stress'] = qm_stress + mm_stress
