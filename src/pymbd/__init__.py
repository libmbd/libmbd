import re
from importlib.metadata import PackageNotFoundError, version

from .pymbd import (
    ang,
    atomic_polarizabilities,
    from_volumes,
    mbd_energy,
    mbd_energy_species,
    molecular_polarizability,
    screening,
    screening_matrix,
)

__all__ = [
    'mbd_energy',
    'mbd_energy_species',
    'screening',
    'screening_matrix',
    'from_volumes',
    'atomic_polarizabilities',
    'molecular_polarizability',
    'ang',
]

try:
    __version__ = version('pymbd')
    __version__ = re.split('[.-]', __version__, maxsplit=3)
    __version__ = (*map(int, __version__[:3]), *__version__[3:])
except PackageNotFoundError:
    __version__ = None
