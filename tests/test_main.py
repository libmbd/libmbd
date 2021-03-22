from pytest import approx


def test_main():
    from pymbd.__main__ import ene, ene_expected

    assert ene == approx(ene_expected, rel=1e-10)
