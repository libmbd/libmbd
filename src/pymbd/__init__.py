import re
from importlib.metadata import PackageNotFoundError, version

from .pymbd import ang, from_volumes, mbd_energy, mbd_energy_species, screening

__all__ = ['mbd_energy', 'mbd_energy_species', 'screening', 'ang', 'from_volumes']
try:
    __version__ = version('pymbd')
    __version__ = re.split('[.-]', __version__, maxsplit=3)
    __version__ = (*map(int, __version__[:3]), *__version__[3:])
except PackageNotFoundError:
    __version__ = None
