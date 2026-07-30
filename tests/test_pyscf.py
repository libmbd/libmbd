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
    # An isolated closed-shell atom has no environment, so its Hirshfeld volume
    # ratio must be exactly 1 -- which holds because the free-atom reference uses
    # the molecular basis.
    mf = dft.RKS(gto.M(atom=f'{elem} 0 0 0', basis='def2-svp', verbose=0), xc='PBE')
    mf.run()
    assert hirshfeld_volume_ratios(mf)[0] == approx(1.0, abs=1e-4)


@pytest.mark.no_scalapack
@pytest.mark.parametrize('basis', ['def2-svp', 'def2-tzvp'])
def test_separated_dimer_ratios_are_one(basis):
    # A well-separated argon dimer is two free atoms, so both ratios must be 1
    # in *any* basis. This fails against a fixed reference basis (0.89 in
    # def2-SVP against aug-cc-pVQZ), which is why the molecular basis is used.
    mf = dft.RKS(gto.M(atom='Ar 0 0 0; Ar 0 0 12.0', basis=basis, verbose=0), xc='PBE')
    mf.run()
    assert hirshfeld_volume_ratios(mf) == approx([1.0, 1.0], abs=1e-4)


@pytest.mark.no_scalapack
def test_water_volume_ratios():
    ratios = hirshfeld_volume_ratios(_water())
    assert ratios == approx([0.95199, 0.549794, 0.549794], abs=1e-3)
    # oxygen stays near-free, hydrogens contract strongly in the O-H bonds
    assert ratios[0] > ratios[1]


@pytest.mark.no_scalapack
def test_water_mbd_energy():
    ene = mbd_energy(_water())
    assert ene == approx(-0.00024035, abs=1e-6)


@pytest.mark.no_scalapack
def test_hartree_fock_is_rejected():
    from pyscf import scf

    mf = scf.RHF(gto.M(atom='Ne 0 0 0', basis='def2-svp', verbose=0)).run()
    with pytest.raises(ValueError, match='KS-DFT'):
        hirshfeld_volume_ratios(mf)


@pytest.mark.no_scalapack
def test_ghost_atoms_are_rejected():
    # counterpoise fragments carry ghost centers, which have no free-atom
    # reference; silently they would give a zero free volume (NaN ratio)
    mol = gto.M(
        atom=[['O', (0, 0, 0)], ['H', (0, 0, 1.8)], ['ghost:O', (0, 0, 5.0)]],
        basis='def2-svp',
        unit='Bohr',
        spin=1,
        verbose=0,
    )
    mf = dft.UKS(mol, xc='PBE').run()
    with pytest.raises(ValueError, match='ghost'):
        hirshfeld_volume_ratios(mf)


@pytest.mark.no_scalapack
def test_mbd_energy_requires_known_beta():
    mf = dft.RKS(gto.M(atom='Ne 0 0 0', basis='def2-svp', verbose=0), xc='M06-L')
    mf.run()
    with pytest.raises(ValueError, match='beta'):
        mbd_energy(mf)
