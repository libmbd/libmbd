import re

import pkg_resources

from .pymbd import (ang, from_volumes, mbd_energy, mbd_energy_species,
                    screening, atomic_polarizability_tensors, molecular_polarizability_tensor)

__all__ = ['mbd_energy', 'mbd_energy_species', 'screening', 'ang',
           'from_volumes', 'atomic_polarizability_tensors', 'molecular_polarizability_tensor']

try:
    __version__ = pkg_resources.get_distribution('pymbd').version
    __version__ = re.split('[.-]', __version__, 3)
    __version__ = (*map(int, __version__[:3]), *__version__[3:])
except pkg_resources.DistributionNotFound:
    __version__ = None
