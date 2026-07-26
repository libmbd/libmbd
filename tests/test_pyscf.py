import pytest
from pytest import approx

pyscf = pytest.importorskip('pyscf')

from pyscf import dft, gto  # noqa: E402

from pymbd.pyscf import hirshfeld_volume_ratios, mbd_energy  # noqa: E402


def _water():
    mol = gto.M(
        atom='O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692',
        basis='def2-svp',
        verbose=0,
    )
    return dft.RKS(mol, xc='PBE').run()


@pytest.mark.no_scalapack
@pytest.mark.parametrize('elem', ['Ne', 'Ar'])
def test_isolated_atom_ratio_is_one(elem):
    # With a matching free-atom basis, an isolated closed-shell atom has no
    # environment and its Hirshfeld volume ratio must be exactly 1.
    mf = dft.RKS(gto.M(atom=f'{elem} 0 0 0', basis='def2-svp', verbose=0), xc='PBE')
    mf.run()
    ratio = hirshfeld_volume_ratios(mf, free_atom_basis='def2-svp')[0]
    assert ratio == approx(1.0, abs=1e-4)


@pytest.mark.no_scalapack
def test_water_volume_ratios():
    ratios = hirshfeld_volume_ratios(_water(), free_atom_basis='def2-svp')
    assert ratios == approx([0.95199, 0.549794, 0.549794], abs=1e-3)
    # oxygen stays near-free, hydrogens contract strongly in the O-H bonds
    assert ratios[0] > ratios[1]


@pytest.mark.no_scalapack
def test_water_mbd_energy():
    ene = mbd_energy(_water(), free_atom_basis='def2-svp')
    assert ene == approx(-0.00024035, abs=1e-6)


@pytest.mark.no_scalapack
def test_mbd_energy_requires_known_beta():
    mf = dft.RKS(gto.M(atom='Ne 0 0 0', basis='def2-svp', verbose=0), xc='M06-L')
    mf.run()
    with pytest.raises(ValueError, match='beta'):
        mbd_energy(mf, free_atom_basis='def2-svp')
