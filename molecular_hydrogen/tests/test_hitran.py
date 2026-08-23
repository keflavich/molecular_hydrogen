"""
Tests that query HITRAN.  They are skipped unless pytest is run with
``--remote-data``.
"""
import numpy as np
import pytest
from astropy import units as u

from molecular_hydrogen import h2
from molecular_hydrogen.tests.test_h2 import KNOWN_WAVELENGTHS

pytestmark = pytest.mark.remote_data


@pytest.fixture(scope='module')
def tbl():
    return h2.get_h2_tbl(min_wavelength=1*u.um, max_wavelength=6*u.um)


class TestTableParsing:

    def test_nonempty(self, tbl):
        assert len(tbl) > 1000

    def test_wavelength_range(self, tbl):
        assert (tbl['wl'].quantity >= 1*u.um).all()
        assert (tbl['wl'].quantity <= 6*u.um).all()

    def test_branches(self, tbl):
        assert set(tbl['branch']) <= set(h2.transdiff)

    def test_transition_types(self, tbl):
        """
        HITRAN flags H2 lines as electric quadrupole (q) or magnetic dipole
        (m); the trailing flag must not be parsed as part of J.
        """
        assert set(tbl['transition_type']) == {'q', 'm'}
        assert (tbl['transition_type'] == 'm').sum() > 0

    def test_magnetic_dipole_lines_parse(self, tbl):
        """
        Regression test: ``int(x.strip('q'))`` raised on quanta like '28m'.
        """
        m = tbl[tbl['transition_type'] == 'm']
        assert len(m) > 0
        assert (m['jl'] >= 0).all()
        # magnetic dipole transitions are Q branch (delta J = 0)
        assert (m['ju'] == m['jl']).all()

    def test_high_j_present(self, tbl):
        """The J values must reach well past 14, not be truncated at 9."""
        assert tbl['jl'].max() > 20

    def test_j_consistency(self, tbl):
        for branch, dj in h2.transdiff.items():
            sel = tbl['branch'] == branch
            if sel.sum():
                assert (tbl['ju'][sel] - tbl['jl'][sel] == dj).all()

    def test_energy_units(self, tbl):
        assert tbl['el_K'].unit == u.K
        assert (np.asarray(tbl['el_K']) >= 0).all()

    def test_wavelength_matches_frequency(self, tbl):
        expected = (tbl['nu']*u.cm**-1).to(u.um, u.spectral())
        np.testing.assert_allclose(tbl['wl'].quantity.value, expected.value)

    def test_linename_column(self, tbl):
        assert '1-0 S(1)' in set(tbl['linename'])

    def test_caching(self):
        first = h2.get_h2_tbl(min_wavelength=2*u.um, max_wavelength=3*u.um)
        second = h2.get_h2_tbl(min_wavelength=2*u.um, max_wavelength=3*u.um)
        assert first is second


class TestWavelengthsAgainstHitran:

    @pytest.mark.parametrize('linename,wavelength',
                             sorted(k for k in KNOWN_WAVELENGTHS.items()
                                    if 1 < k[1] < 6))
    def test_dabrowski_matches_hitran(self, tbl, linename, wavelength):
        """
        The rest wavelengths computed from the Dabrowski (1984) term values
        must agree with HITRAN.  ``extended=False`` is essential: the extended
        table is itself derived from HITRAN, so the comparison would otherwise
        be circular.
        """
        row = tbl[tbl['linename'] == linename]
        assert len(row) > 0, "%s missing from HITRAN table" % linename
        ours = h2.linename_to_restwl(linename, extended=False)
        np.testing.assert_allclose(ours.to(u.um).value, float(row['wl'][0]),
                                   rtol=2e-4)

    def test_all_tabulated_levels_agree_with_hitran(self, tbl):
        """
        Cross-check every transition we can compute from `resten` against
        HITRAN.  This is what caught the V=0, J=12 typo.  ``extended=False``
        keeps this a real check on the Dabrowski table rather than a
        comparison of HITRAN with itself.
        """
        quad = tbl[tbl['transition_type'] == 'q']
        checked = 0
        for row in quad:
            vu, vl = row['global_upper_quanta'], row['global_lower_quanta']
            ju, jl = row['ju'], row['jl']
            if (vu >= h2.resten.shape[1] or vl >= h2.resten.shape[1] or
                    ju >= h2.resten.shape[0] or jl >= h2.resten.shape[0] or
                    ju < 0 or jl < 0):
                continue
            ours = h2.restwl(vu, vl, ju, jl,
                             extended=False).to(u.um).value
            theirs = float(row['wl'])
            assert abs(ours - theirs) / theirs < 5e-4, \
                "%s: %.6f vs HITRAN %.6f" % (row['linename'], ours, theirs)
            checked += 1
        assert checked > 100


class TestAval:

    @pytest.mark.parametrize('linename,expected', [
        ('1-0 S(1)', 3.47e-7),
        ('1-0 S(0)', 2.53e-7),
        ('1-0 S(5)', 3.96e-7),
        ('1-0 Q(1)', 4.29e-7),
        ('1-0 Q(3)', 2.78e-7),
        ('2-1 S(1)', 4.98e-7),
        ('3-2 S(1)', 5.14e-7),
    ])
    def test_against_ukirt_table(self, linename, expected):
        """HITRAN A values must match the published UKIRT list to ~1%."""
        vu, vl = linename.split()[0].split('-')
        branch = linename.split()[1][0]
        jl = int(linename.split('(')[1].rstrip(')'))
        ju = jl + h2.transdiff[branch]
        A = h2.aval(int(vu), int(vl), ju, jl)
        assert A.unit == u.s**-1
        np.testing.assert_allclose(A.value, expected, rtol=0.01)

    def test_quadrupole_is_the_default(self):
        """
        Regression test: the Q-branch lines are listed twice in HITRAN, once
        as quadrupole and once as magnetic dipole.  Taking the first match
        returned the magnetic dipole rate for 1-0 Q(3), 100x too small.
        """
        quad = h2.aval(1, 0, 3, 3)
        mag = h2.aval(1, 0, 3, 3, transition_type='m')
        np.testing.assert_allclose(quad.value, 2.78e-7, rtol=0.01)
        assert mag < quad / 10

    def test_missing_transition_raises(self):
        with pytest.raises(ValueError):
            h2.aval(1, 0, 3, 1, transition_type='m')  # S branch has no m lines


class TestEmission:

    def test_emission_per_atom_oneline(self):
        n = h2.emission_per_atom_oneline(2000., 3., 1, 3, 1)
        assert n.unit == u.s**-1
        assert n.value > 0

    def test_emission_rises_with_temperature(self):
        cold = h2.emission_per_atom_oneline(500., 3., 1, 3, 1)
        hot = h2.emission_per_atom_oneline(3000., 3., 1, 3, 1)
        assert hot > cold

    def test_vectorized(self):
        temps = np.array([1000., 2000., 3000.])
        result = h2.vectorized_emission_per_atom_oneline(temps, 3., 1, 3, 1)
        assert result.shape == (3,)
        assert (np.diff(result.value) > 0).all()

    def test_emission_per_atom(self):
        n = h2.emission_per_atom(2000., 3.)
        assert len(n) > 0
        assert all(v.value >= 0 for v in n.values())


class TestExtendedTableAgainstHitran:

    def test_extended_wavelengths_reproduce_hitran(self, tbl):
        """
        The bundled level energies must reproduce the HITRAN line wavelengths
        they were derived from, for every transition in the table.
        """
        quad = tbl[tbl['transition_type'] == 'q']
        energies = h2.get_level_energies()
        checked = 0
        for row in quad:
            vu, vl = row['global_upper_quanta'], row['global_lower_quanta']
            ju, jl = row['ju'], row['jl']
            if (ju < 0 or jl < 0 or ju >= energies.shape[0] or
                    jl >= energies.shape[0] or vu >= energies.shape[1] or
                    vl >= energies.shape[1]):
                continue
            if not (np.isfinite(energies[ju][vu]) and
                    np.isfinite(energies[jl][vl])):
                continue
            ours = h2.restwl(vu, vl, ju, jl).to(u.um).value
            theirs = float(row['wl'])
            assert abs(ours - theirs) / theirs < 1e-5, \
                "%s: %.6f vs HITRAN %.6f" % (row['linename'], ours, theirs)
            checked += 1
        assert checked > 1000

    def test_every_hitran_level_is_in_the_table(self, tbl):
        """No level appearing in the line list may be missing from the grid."""
        energies = h2.get_level_energies()
        for row in tbl:
            for V, J in ((row['global_lower_quanta'], row['jl']),
                         (row['global_upper_quanta'], row['ju'])):
                assert V < energies.shape[1] and J < energies.shape[0], \
                    "level V=%s J=%s outside the grid" % (V, J)
                assert np.isfinite(energies[J][V]), \
                    "level V=%s J=%s missing from the grid" % (V, J)
