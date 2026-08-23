"""
Molecular hydrogen rovibrational line properties.

Rest wavelengths come from the measured X-state term values of
Dabrowski (1984), Can. J. Phys. 62, 1639.  Einstein A coefficients are looked
up from HITRAN (molecule 45) via `astroquery`.
"""
import re

import numpy as np
from astropy import units as u
from astropy import constants

# define physical constants to high precision
h = constants.h.cgs.value  # 6.626068e-27
c = constants.c.cgs.value  # 2.99792e10
kb = constants.k_B.cgs.value  # 1.3806503e-16
e = constants.e.esu  # 4.803e-12

# Dunham-like coefficients for the X state, from NIST
We = 4401.21
Be = 60.853
WeXe = 121.33
De = .0471
Ae = 3.062

# http://www.nist.gov/data/nsrds/NSRDS-NBS31.pdf
dissociation_energy = (432*u.kJ / constants.N_A.value).to(u.eV)
dissociation_temperature = (dissociation_energy/constants.k_B).to(u.K)


def pyh2level_energy(V, J):
    """ Returns the theoretical level energy as a function of the
    vibrational (V) and rotational (J) state of the molecule.

    Constants are from NIST:
    http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1000#Diatomic
    (see the bottom of the table)

    Returns a value in ergs

    (8x slower than cython version)

    Notes
    -----
    This is a low-order expansion that is only reliable for low V and J; for
    quantitative work use the measured term values in `resten` (see
    `level_energy`) instead.
    """

    # not used re=0.74144

    return (h * c *
            (We*(V+0.5) + Be*(J*(J+1)) - WeXe*(V+.5)**2 - De*J**2*(J+1)**2 - Ae*(V+.5)*(J+1)*J))


def pylevel_population(temperature, orthopararatio, jmax, vmax):
    """
    Pure-python version of `level_population`.

    Compute LTE level populations, normalized to sum to one over the
    ``0 <= V < vmax``, ``0 <= J < jmax`` levels included.
    """

    para = 1./(1.+orthopararatio)
    ortho = 1. - para

    kbt = kb * temperature
    if kbt == 0:
        # can't have div by zero...
        return {(V, J): np.nan
                for J in range(0, jmax)
                for V in range(0, vmax)}

    # Compute thermal level populations
    lp = {(V, J):
          (2*J+1)*np.exp(-h2level_energy(V, J)/kbt) *
          (para if J % 2 == 0 else ortho)
          for J in range(0, jmax)
          for V in range(0, vmax)}
    total_population = sum(lp.values())
    if total_population == 0:
        return lp
    return {k: v/total_population for k, v in lp.items()}


# The cython implementations are ~8x faster, but they require a working
# compiler.  Fall back to the pure-python versions when that is unavailable so
# that importing the package never hard-fails.
try:
    import pyximport
    pyximport.install()
    from .h2level_energy import h2level_energy
    from .level_population import level_population
    HAS_COMPILED_EXTENSIONS = True
except Exception:  # pragma: no cover - depends on the build environment
    h2level_energy = pyh2level_energy
    level_population = pylevel_population
    HAS_COMPILED_EXTENSIONS = False


# Measured X-state term values in cm^-1, indexed as resten[J][V].
# Dabrowski, I. 1984, Can. J. Phys. 62, 1639 ("The Lyman and Werner bands of
# H2"), Table 3.  Values agree with the HITRAN level energies to <0.2 cm^-1.
resten = np.array(
      [[   0.00,   4161.14,   8086.93,  11782.36,  15250.31,  18491.92],
       [ 118.50,   4273.75,   8193.81,  11883.51,  15345.81,  18581.71],
       [ 354.35,   4497.82,   8406.29,  12084.66,  15535.7 ,  18760.28],
       [ 705.54,   4831.41,   8722.7 ,  12384.14,  15818.27,  19026.01],
       [1168.78,   5271.36,   9139.86,  12778.78,  16190.66,  19375.99],
       [1740.21,   5813.95,   9654.15,  13265.27,  16649.48,  19807.03],
       [2414.76,   6454.28,  10261.2 ,  13839.18,  17190.36,  20314.77],
       [3187.57,   7187.44,  10955.68,  14495.46,  17808.76,  20894.94],
       [4051.73,   8007.77,  11732.12,  15228.88,  18499.06,  21542.14],
       [5001.97,   8908.28,  12584.8 ,  16033.83,  19256.43,  22251.21],
       [6030.81,   9883.79,  13507.42,  16904.02,  20074.45,  23016.21],
       [7132.03,  10927.12,  14493.58,  17833.83,  20947.48,  23831.68],
       [8298.60,  12031.44,  15537.15,  18816.78,  21869.48,  24691.77],
       [9523.82,  13191.06,  16632.1 ,  19847.08,  22834.61,  25590.22]])


def level_energy(V, J):
    """
    Measured term value of the X-state (V, J) level.

    Parameters
    ----------
    V, J : int
        Vibrational and rotational quantum numbers

    Returns
    -------
    Energy above the (V=0, J=0) level, in cm^-1
    """
    if not (0 <= J < resten.shape[0]) or not (0 <= V < resten.shape[1]):
        raise NotImplementedError("No data for V=%s, J=%s" % (V, J))
    return resten[J][V] * u.cm**-1


def restwl(vu, vl, ju, jl):
    """
    Rest wavelength of a rovibrational transition.

    Uses energy levels measured by Dabrowski, Can. J. Phys. 62, 1639, 1984.

    Parameters
    ----------
    vu,vl : int,int
        upper and lower vibrational states
    ju,jl : int,int
        upper and lower rotational states

    Returns
    -------
    Wavelength in microns

    Examples
    --------
    >>> round(float(restwl(1, 0, 3, 1).value), 6)
    2.121831
    """
    for V, J, which in ((vu, ju, 'upper'), (vl, jl, 'lower')):
        if not (0 <= J < resten.shape[0]) or not (0 <= V < resten.shape[1]):
            raise NotImplementedError("No data for the %s level v=%s, J=%s"
                                      % (which, V, J))

    denergy = resten[ju][vu] - resten[jl][vl]
    if denergy <= 0:
        raise ValueError("The upper level (v=%s, J=%s) is not above the lower "
                         "level (v=%s, J=%s)" % (vu, ju, vl, jl))

    # resten is in inverse cm.  Convert to m
    dl = .01/denergy
    return (dl * u.m).to(u.um)


transdiff = {'S': 2, 'R': 1, 'Q': 0, 'P': -1, 'O': -2}
reverse_transdiff = {v: k for k, v in transdiff.items()}

_linename_re = re.compile(r'([0-9]+)\s*-\s*([0-9]+)')
_branch_re = re.compile(r'([%s])\(([0-9]+)\)' % ''.join(sorted(transdiff)))


def linename_to_restwl(linename):
    """
    Parse a line name of the form S(1) 1-0, Q(2) 2-1, or 0-0 S(0), etc.

    Parameters
    ----------
    linename : str
        String of the form %i-%i %s(%i), where the first numbers are the
        vibrational levels, the letter is the rotational branch (one of
        O, P, Q, R, S), and the last number is the lower rotational state of
        the transition

    Returns
    -------
    Wavelength in microns

    Examples
    --------
    >>> round(float(linename_to_restwl('1-0 S(1)').value), 6)
    2.121831
    >>> round(float(linename_to_restwl('1-0 Q(1)').value), 6)
    2.406594
    """
    vmatch = _linename_re.search(linename)
    bmatch = _branch_re.search(linename)
    if vmatch is None or bmatch is None:
        raise ValueError("Could not parse line name %r; expected something "
                         "like '1-0 S(1)'" % (linename,))
    upper, lower = vmatch.groups()
    transtype, jl = bmatch.groups()
    ju = int(jl) + transdiff[transtype]
    return restwl(int(upper), int(lower), ju, int(jl))


def restwl_to_linename(vu, vl, ju, jl):
    """
    Inverse of the name parsing in `linename_to_restwl`: build the canonical
    line name for a transition.

    Examples
    --------
    >>> restwl_to_linename(1, 0, 3, 1)
    '1-0 S(1)'
    """
    if ju - jl not in reverse_transdiff:
        raise ValueError("No branch name for ju-jl=%s" % (ju - jl))
    return '%i-%i %s(%i)' % (vu, vl, reverse_transdiff[ju - jl], jl)


# Einstein A values from
# http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
# retained for reference and for testing the HITRAN lookup; `aval` is the
# preferred interface.
aval_dict = {1: {'S0': 2.53e-7,
                 'S1': 3.47e-7,
                 'S2': 3.98e-7,
                 'S3': 4.21e-7,
                 'S4': 4.19e-7,
                 'S5': 3.96e-7,
                 'S6': 3.54e-7,
                 'S7': 2.98e-7,
                 'S8': 2.34e-7,
                 'S9': 1.68e-7,
                 'Q1': 4.29e-7,
                 'Q2': 3.03e-7,
                 'Q3': 2.78e-7,
                 'Q4': 2.65e-7},
             2: {'S0': 3.68e-7,
                 'S1': 4.98e-7,
                 'S2': 5.60e-7,
                 'S3': 5.77e-7,
                 'S4': 5.57e-7, },
             3: {'S0': 3.88e-7,
                 'S1': 5.14e-7,
                 'S2': 5.63e-7,
                 'S3': 5.63e-7,
                 'S4': 5.22e-7, }
             }

# HITRAN queries are slow; keep the retrieved tables around.
_h2_tbl_cache = {}


def get_h2_tbl(min_wavelength=0.9*u.um, max_wavelength=30*u.um, cache=True):
    """
    Retrieve the HITRAN H2 line list over a wavelength range.

    Parameters
    ----------
    min_wavelength, max_wavelength : `~astropy.units.Quantity`
        Wavelength limits of the query
    cache : bool
        Re-use a previously retrieved table for the same wavelength range?

    Returns
    -------
    `~astropy.table.Table` with the HITRAN columns plus the decoded quantum
    numbers ``branch``, ``jl``, ``ju``, ``transition_type``, the wavelength
    ``wl``, the lower level energy ``el_K``, and the ``linename``.

    Notes
    -----
    The HITRAN local quanta strings look like ``'     S  1q     '``: the letter
    is the rotational branch, the number is the lower rotational state, and the
    trailing flag is the transition type (``q`` = electric quadrupole,
    ``m`` = magnetic dipole).  Both are tabulated for the Q branch, and the
    magnetic dipole rates are ~2-3 orders of magnitude smaller.
    """
    key = (min_wavelength.to(u.um).value, max_wavelength.to(u.um).value)
    if cache and key in _h2_tbl_cache:
        return _h2_tbl_cache[key]

    from astroquery.hitran import Hitran
    tbl = Hitran.query_lines(molecule_number=45,
                             isotopologue_number=1,
                             min_frequency=max_wavelength.to(u.THz, u.spectral()),
                             max_frequency=min_wavelength.to(u.THz, u.spectral()))
    tbl['branch'] = [x[5] for x in tbl['local_lower_quanta']]
    # the trailing character is the transition type (q or m), not part of J
    tbl['jl'] = [int(x.split()[1].rstrip('qm')) for x in tbl['local_lower_quanta']]
    tbl['transition_type'] = [x.split()[1][-1] for x in tbl['local_lower_quanta']]
    tbl['ju'] = [LL + transdiff[BB] for LL, BB in zip(tbl['jl'], tbl['branch'])]
    tbl['wl'] = (tbl['nu']*u.cm**-1).to(u.um, u.spectral())
    tbl['el_K'] = (tbl['elower']*u.cm**-1 * constants.c * constants.h / constants.k_B).to(u.K)
    tbl['global_lower_quanta'] = tbl['global_lower_quanta'].astype('int')
    tbl['global_upper_quanta'] = tbl['global_upper_quanta'].astype('int')
    tbl['linename'] = [restwl_to_linename(VU, VL, JU, JL)
                       for VU, VL, JU, JL in zip(tbl['global_upper_quanta'],
                                                 tbl['global_lower_quanta'],
                                                 tbl['ju'], tbl['jl'])]

    if cache:
        _h2_tbl_cache[key] = tbl
    return tbl


def aval(vu, vl, ju, jl, transition_type='q'):
    """
    Einstein A coefficient of a rovibrational transition, from HITRAN.

    Parameters
    ----------
    vu,vl : int,int
        upper and lower vibrational states
    ju,jl : int,int
        upper and lower rotational states
    transition_type : 'q' or 'm'
        HITRAN tabulates both electric quadrupole ('q') and magnetic dipole
        ('m') rates for the Q branch.  The quadrupole lines are the ones of
        astrophysical interest and are the default; the magnetic dipole rates
        are 2-3 orders of magnitude smaller.

    Returns
    -------
    A in s^-1
    """
    tbl = get_h2_tbl()
    match = ((tbl['jl'] == jl) & (tbl['ju'] == ju) &
             (tbl['global_lower_quanta'] == vl) &
             (tbl['global_upper_quanta'] == vu) &
             (tbl['transition_type'] == transition_type))
    if not match.any():
        raise ValueError("Transition %s is not in the HITRAN table"
                         % restwl_to_linename(vu, vl, ju, jl))
    return tbl['a'][match][0] * u.s**-1


def emission_per_atom(temperature, orthopararatio, jmax=14, vmax=6):
    """
    Compute the emission per H2 molecule in the rovibrational states assuming
    thermal level populations given a temperature and an ortho/para ratio of H2

    Returns
    -------
    dict mapping (V, J) of the upper level of each S-branch transition to the
    number of photons emitted per molecule per second
    """
    lp = level_population(temperature, orthopararatio, jmax, vmax)

    nphotons = {}
    for V in range(1, vmax):
        for J in range(2, jmax):
            try:
                A = aval(V, V-1, J, J-2)
            except ValueError:
                continue
            nphotons[(V, J)] = A * lp[(V, J)]

    return nphotons


def emission_per_atom_oneline(temperature, orthopararatio, vu, ju, jl,
                              vl=None, jmax=14, vmax=6):
    """
    Compute the emission per H2 molecule in a single rovibrational line
    assuming thermal level populations given a temperature and an ortho/para
    ratio of H2

    Parameters
    ----------
    vl : int or None
        Lower vibrational state.  Defaults to ``vu - 1``.
    """
    if vl is None:
        vl = vu - 1

    lp = level_population(temperature, orthopararatio, jmax, vmax)

    return aval(vu, vl, ju, jl) * lp[(vu, ju)]


def vectorized_emission_per_atom_oneline(temperatures, orthopararatio, vu, ju,
                                         jl, vl=None, jmax=14, vmax=6):
    return u.Quantity([emission_per_atom_oneline(t, orthopararatio, vu, ju, jl,
                                                 vl=vl, jmax=jmax, vmax=vmax)
                       for t in temperatures])
