import numpy as np
import pytest
from astropy import units as u

from molecular_hydrogen import h2


# Rest wavelengths in vacuum, in microns, from HITRAN (molecule 45).
# `restwl` reproduces these from the Dabrowski (1984) term values; the term
# values are quoted to 0.01 cm^-1, which limits the agreement to ~1e-4
# fractionally (worst for the low-energy pure-rotational lines).
KNOWN_WAVELENGTHS = {
    '0-0 S(0)': 28.21884,
    '0-0 S(1)': 17.03485,
    '0-0 S(2)': 12.27861,
    '0-0 S(3)': 9.66491,
    '0-0 S(5)': 6.90951,
    '1-0 S(0)': 2.22329,
    '1-0 S(1)': 2.12183,
    '1-0 S(2)': 2.03376,
    '1-0 S(3)': 1.95756,
    '1-0 Q(1)': 2.40659,
    '1-0 Q(2)': 2.41344,
    '1-0 Q(3)': 2.42373,
    '1-0 O(2)': 2.62688,
    '1-0 O(3)': 2.80252,
    '2-1 S(1)': 2.24772,
    '2-1 S(2)': 2.15422,
    '3-2 S(3)': 2.20140,
    '2-0 S(0)': 1.18958,
    '2-0 S(1)': 1.16222,
}


class TestRestWavelengths:

    @pytest.mark.parametrize('linename,wavelength', sorted(KNOWN_WAVELENGTHS.items()))
    def test_known_lines(self, linename, wavelength):
        result = h2.linename_to_restwl(linename)
        assert result.unit.is_equivalent(u.um)
        # limited by the 0.01 cm^-1 precision of the tabulated term values
        np.testing.assert_allclose(result.to(u.um).value, wavelength, rtol=2e-4)

    def test_returns_quantity(self):
        assert h2.restwl(1, 0, 3, 1).unit == u.um

    def test_name_order_does_not_matter(self):
        assert (h2.linename_to_restwl('1-0 S(1)') ==
                h2.linename_to_restwl('S(1) 1-0'))

    @pytest.mark.parametrize('branch,dj', [('O', -2), ('P', -1), ('Q', 0),
                                           ('R', 1), ('S', 2)])
    def test_all_branches_parse(self, branch, dj):
        """Every branch in ``transdiff`` must be accepted by the parser."""
        jl = 3
        expected = h2.restwl(1, 0, jl + dj, jl)
        assert h2.linename_to_restwl('1-0 %s(%i)' % (branch, jl)) == expected

    @pytest.mark.parametrize('bad', ['', 'S(1)', '1-0', 'gibberish',
                                     '1-0 Z(1)', '1-0 S1'])
    def test_unparseable_names_raise(self, bad):
        with pytest.raises(ValueError):
            h2.linename_to_restwl(bad)

    def test_out_of_range_raises(self):
        with pytest.raises(NotImplementedError):
            h2.restwl(1, 0, 40, 38)
        with pytest.raises(NotImplementedError):
            h2.restwl(20, 0, 3, 1)
        # a corner of the grid that HITRAN does not cover
        with pytest.raises(NotImplementedError):
            h2.restwl(14, 0, 31, 29)

    def test_out_of_range_raises_without_extended_table(self):
        with pytest.raises(NotImplementedError):
            h2.restwl(1, 0, 30, 28, extended=False)
        with pytest.raises(NotImplementedError):
            h2.restwl(9, 0, 3, 1, extended=False)

    def test_extended_table_is_used_by_default(self):
        """
        The extended (HITRAN) term values are quoted to 0.0001 cm^-1 against
        0.01 cm^-1 for `resten`, so the two differ slightly.
        """
        assert (h2.restwl(1, 0, 3, 1) !=
                h2.restwl(1, 0, 3, 1, extended=False))
        np.testing.assert_allclose(h2.restwl(1, 0, 3, 1).value,
                                   h2.restwl(1, 0, 3, 1, extended=False).value,
                                   rtol=1e-5)

    def test_negative_quantum_number_raises(self):
        # an O-branch line out of the bottom of the ladder
        with pytest.raises(NotImplementedError):
            h2.linename_to_restwl('1-0 O(1)')

    def test_non_emitting_transition_raises(self):
        """The 'upper' level must actually lie above the lower level."""
        with pytest.raises(ValueError):
            h2.restwl(0, 1, 1, 3)
        with pytest.raises(ValueError):
            h2.restwl(0, 0, 1, 1)

    def test_roundtrip_name(self):
        for vu, vl, ju, jl in [(1, 0, 3, 1), (2, 1, 2, 0), (1, 0, 1, 3)]:
            name = h2.restwl_to_linename(vu, vl, ju, jl)
            assert h2.linename_to_restwl(name) == h2.restwl(vu, vl, ju, jl)


class TestLevelEnergies:

    def test_shape(self):
        assert h2.resten.shape == (14, 6)

    def test_ground_state_is_zero(self):
        assert h2.resten[0][0] == 0

    def test_monotonic_in_J(self):
        assert (np.diff(h2.resten, axis=0) > 0).all()

    def test_monotonic_in_V(self):
        assert (np.diff(h2.resten, axis=1) > 0).all()

    def test_rotational_ladder_is_smooth(self):
        """
        Coarse sanity check on the shape of each rotational ladder: successive
        second differences of E(J) vary by <10% for the real term values, so a
        20% threshold flags gross corruption without being brittle.
        """
        for V in range(h2.resten.shape[1]):
            d2 = np.diff(h2.resten[:, V], n=2)
            assert (np.abs(np.diff(d2)) < 0.2 * np.abs(d2[:-1])).all(), \
                "kink in the V=%i rotational ladder: %s" % (V, d2)

    def test_v0_j12_value(self):
        """
        Regression test: this entry was 8296.61, ~2 cm^-1 below the HITRAN and
        Dabrowski (1984) value, which shifted 0-0 S(10) by ~4 nm.
        """
        np.testing.assert_allclose(h2.resten[12][0], 8298.60, atol=0.2)

    def test_level_energy_units(self):
        assert h2.level_energy(0, 1).unit == u.cm**-1
        np.testing.assert_allclose(h2.level_energy(0, 1).value, 118.4869)
        np.testing.assert_allclose(
            h2.level_energy(0, 1, extended=False).value, 118.50)

    def test_level_energy_out_of_range(self):
        with pytest.raises(NotImplementedError):
            h2.level_energy(0, 40)
        with pytest.raises(NotImplementedError):
            h2.level_energy(0, 20, extended=False)


class TestLevelPopulation:

    def test_normalized(self):
        lp = h2.level_population(2000., 3., 14, 6)
        np.testing.assert_allclose(sum(lp.values()), 1.0)

    def test_keys(self):
        lp = h2.level_population(2000., 3., 5, 3)
        assert set(lp) == {(V, J) for V in range(3) for J in range(5)}

    def test_zero_temperature_is_nan(self):
        lp = h2.level_population(0., 3., 5, 3)
        assert all(np.isnan(v) for v in lp.values())

    def test_ortho_para_ratio(self):
        """
        Summing over J of each parity must recover the requested OPR (the
        even-J/para and odd-J/ortho weights are applied before normalization).
        """
        opr = 3.
        lp = h2.level_population(5000., opr, 14, 6)
        para = sum(v for (V, J), v in lp.items() if J % 2 == 0)
        ortho = sum(v for (V, J), v in lp.items() if J % 2 == 1)
        # not exactly opr, since the partition functions differ, but the
        # weights must at least put more molecules in ortho for opr > 1
        assert ortho > para

    def test_hotter_populates_higher_levels(self):
        cold = h2.level_population(500., 3., 14, 6)
        hot = h2.level_population(5000., 3., 14, 6)
        assert hot[(1, 3)] > cold[(1, 3)]

    def test_cython_matches_python(self):
        cy = h2.level_population(2000., 3., 14, 6)
        py = h2.pylevel_population(2000., 3., 14, 6)
        assert set(cy) == set(py)
        for key in cy:
            np.testing.assert_allclose(cy[key], py[key], rtol=1e-10)


class TestLevelEnergyFormula:

    def test_cython_matches_python(self):
        for V in range(6):
            for J in range(14):
                np.testing.assert_allclose(h2.h2level_energy(V, J),
                                           h2.pyh2level_energy(V, J),
                                           rtol=1e-10)

    def test_approximates_measured_energies(self):
        """
        The Dunham expansion is only approximate, but the (V,J) -> (0,0)
        differences should track the measured term values to a few percent for
        the low-lying levels it is fit to.
        """
        e00 = h2.pyh2level_energy(0, 0)
        for V in range(3):
            for J in range(1, 8):
                approx = ((h2.pyh2level_energy(V, J) - e00) * u.erg).to(
                    u.cm**-1, u.spectral())
                measured = h2.level_energy(V, J)
                np.testing.assert_allclose(approx.value, measured.value,
                                           rtol=0.05)


class TestParsingHelpers:

    def test_transdiff_covers_all_branches(self):
        assert h2.transdiff == {'S': 2, 'R': 1, 'Q': 0, 'P': -1, 'O': -2}

    def test_reverse_transdiff(self):
        for branch, dj in h2.transdiff.items():
            assert h2.reverse_transdiff[dj] == branch

    def test_linename_construction(self):
        assert h2.restwl_to_linename(1, 0, 3, 1) == '1-0 S(1)'
        assert h2.restwl_to_linename(2, 1, 0, 2) == '2-1 O(2)'
        assert h2.restwl_to_linename(1, 0, 1, 1) == '1-0 Q(1)'

    def test_linename_bad_branch(self):
        with pytest.raises(ValueError):
            h2.restwl_to_linename(1, 0, 6, 1)


def test_dissociation_energy():
    """H2 dissociates at 4.48 eV / 52000 K."""
    np.testing.assert_allclose(h2.dissociation_energy.to(u.eV).value, 4.48,
                               rtol=0.01)
    np.testing.assert_allclose(h2.dissociation_temperature.to(u.K).value,
                               52000, rtol=0.01)


class TestExtendedLevelEnergies:
    """
    The Dabrowski (1984) table stops at V=5, J=13.  `get_level_energies`
    extends this to V=0-14, J=0-31 using the level energies implied by the
    HITRAN line list.
    """

    def test_shape(self):
        assert h2.get_level_energies().shape == (32, 15)

    def test_cached(self):
        assert h2.get_level_energies() is h2.get_level_energies()

    def test_ground_state_is_zero(self):
        assert h2.get_level_energies()[0][0] == 0

    def test_v0_covers_j_up_to_31(self):
        """The specific gap this table was added to fill: v=0, J>13."""
        v0 = h2.get_level_energies()[:, 0]
        assert np.isfinite(v0).all()
        assert len(v0) == 32

    def test_ragged_grid_has_nans(self):
        """High J is only reached for low V; those entries must be NaN."""
        energies = h2.get_level_energies()
        assert np.isnan(energies[31][1])
        assert np.isnan(energies[20][10])

    def test_monotonic_in_J_where_defined(self):
        energies = h2.get_level_energies()
        for V in range(energies.shape[1]):
            column = energies[:, V]
            column = column[np.isfinite(column)]
            assert (np.diff(column) > 0).all(), "V=%i not monotonic" % V

    def test_monotonic_in_V_where_defined(self):
        energies = h2.get_level_energies()
        for J in range(energies.shape[0]):
            row = energies[J, :]
            row = row[np.isfinite(row)]
            assert (np.diff(row) > 0).all(), "J=%i not monotonic" % J

    def test_below_dissociation(self):
        """
        Every tabulated level must lie below the dissociation limit, and the
        topmost one (V=14) should sit just under it: D0 = 36118.07 cm^-1
        (Liu et al. 2009, J. Chem. Phys. 130, 174306).  ``h2.dissociation_energy``
        is derived from a rounded 432 kJ/mol and is ~6 cm^-1 low, so the
        accurate value is used here.
        """
        D0 = 36118.07
        energies = h2.get_level_energies()
        finite = energies[np.isfinite(energies)]
        assert finite.max() <= D0
        assert finite.max() > 0.99 * D0
        np.testing.assert_allclose(
            h2.dissociation_energy.to(u.cm**-1, u.spectral()).value, D0,
            rtol=1e-3)

    def test_agrees_with_dabrowski_table(self):
        """
        Where the two overlap they must agree; `resten` is quoted to 0.01
        cm^-1.
        """
        energies = h2.get_level_energies()
        for J in range(h2.resten.shape[0]):
            for V in range(h2.resten.shape[1]):
                np.testing.assert_allclose(
                    energies[J][V], h2.resten[J][V], atol=0.25,
                    err_msg="V=%i, J=%i" % (V, J))

    def test_high_j_pure_rotational_lines(self):
        """
        Pure rotational 0-0 S(J) lines for J > 13, which could not be computed
        before.  Values from HITRAN.
        """
        for linename, wavelength in [('0-0 S(14)', 3.7245),
                                     ('0-0 S(15)', 3.6263),
                                     ('0-0 S(20)', 3.3812),
                                     ('0-0 S(25)', 3.4111),
                                     ('0-0 S(29)', 3.5992)]:
            result = h2.linename_to_restwl(linename)
            np.testing.assert_allclose(result.to(u.um).value, wavelength,
                                       rtol=1e-4, err_msg=linename)

    def test_rotational_ladder_turns_over(self):
        """
        The 0-0 S(J) wavelengths shorten with J up to about J=20 and then
        lengthen again as the centrifugal distortion flattens the ladder.
        """
        wl = np.array([h2.linename_to_restwl('0-0 S(%i)' % J).value
                       for J in range(10, 30)])
        assert wl.argmin() not in (0, len(wl)-1)
