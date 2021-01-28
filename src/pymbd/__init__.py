import pkg_resources

from .pymbd import ang, from_volumes, mbd_energy, mbd_energy_species, screening

__all__ = ['mbd_energy', 'mbd_energy_species', 'screening', 'ang', 'from_volumes']
try:
    __version__ = tuple(
        map(int, pkg_resources.get_distribution('pymbd').version.split('.')[:3])
    )
except pkg_resources.DistributionNotFound:
    __version__ = None
