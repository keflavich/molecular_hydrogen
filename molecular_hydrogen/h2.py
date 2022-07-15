import numpy as np
import re
from astropy import units as u
from astropy import constants

# define physical constants to high precision
h=constants.h.cgs.value #6.626068e-27
c=constants.c.cgs.value #2.99792e10
k=constants.k_B.cgs # 1.3806503e-16
kb = constants.k_B.cgs.value
e=constants.e.esu #4.803e-12

We=4401.21
Be=60.853
WeXe=121.33
De=.0471
Ae=3.062

# http://www.nist.gov/data/nsrds/NSRDS-NBS31.pdf
dissociation_energy = (432*u.kJ / constants.N_A.value).to(u.eV)
dissociation_temperature = (dissociation_energy/constants.k_B).to(u.K)

import pyximport
pyximport.install()
from .h2level_energy import h2level_energy
from .level_population import level_population

def pyh2level_energy(V,J):
    """ Returns the theoretical level energy as a function of the
    vibrational (V) and rotational (J) state of the molecule.
    
    Constants are from NIST:
    http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1000#Diatomic
    (see the bottom of the table)

    Returns a value in ergs

    (8x slower than cython version)
    """

    # not used re=0.74144

    return (h * c *
            (We*(V+0.5) + Be*(J*(J+1)) - WeXe*(V+.5)**2 - De*J**2*(J+1)**2 - Ae*(V+.5)*(J+1)*J))

# Dalgarno 1984 table
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
       [8296.61,  12031.44,  15537.15,  18816.78,  21869.48,  24691.77],
       [9523.82,  13191.06,  16632.1 ,  19847.08,  22834.61,  25590.22]])


def restwl(vu,vl,ju,jl):
    """
    Uses energy levels measured by Dalgarno, Can J. Physics, 62,1639,1984

    Parameters
    ----------
    vu,vl : int,int
        upper and lower vibrational states
    ju,jl : int,int
        upper and lower rotational states
    Returns
    -------
    Wavelength in microns
    """
    if ju >= resten.shape[0] or vu >= resten.shape[1]:
        raise NotImplementedError("No data for ju=%i,vu=%i" % (ju,vu))
    # resten is in inverse cm.  Convert to m
    dl = .01/(resten[ju][vu]-resten[jl][vl])
    return (dl * u.m).to(u.um)

transdiff = {'S':2,'Q':0,'O':-2}
transdiff = {'S': 2, 'R': 1, 'Q': 0, 'P': -1, 'O': -2}

def linename_to_restwl(linename):
    """
    Parse a line name of the form S(1) 1-0, Q(2) 2-1, or 0-0 S(0), etc.

    Parameters
    ----------
    linename : str
        String of the form %i-%i %s(%i), where the first numbers are the
        vibrational levels, the letter is S, Q, or O indicating the rotational
        branch, and the last number is the lower rotational state of the
        transition
    
    Returns
    -------
    Wavelength in microns

    Examples
    --------
    >>> linename_to_restwl('1-0 S(1)') 
    2.1218313101671793
    >>> linename_to_restwl('1-0 Q(1)')
    2.406594067745623
    """
    upper,lower = re.compile('([0-9]*)-([0-9]*)').search(linename).groups()
    transtype,jl = re.compile('([SQO])\(([0-9]*)\)').search(linename).groups()
    ju = int(jl) + transdiff[transtype]
    rwl = restwl(int(upper),int(lower),ju,int(jl))
    return rwl


aval_dict = {1: {'S0':2.53e-7,
                 'S1':3.47e-7,
                 'S2':3.98e-7,
                 'S3':4.21e-7,
                 'S4':4.19e-7,
                 'S5':3.96e-7,
                 'S6':3.54e-7,
                 'S7':2.98e-7,
                 'S8':2.34e-7,
                 'S9':1.68e-7,
                 'Q1':4.29e-7,
                 'Q2':3.03e-7,
                 'Q3':2.78e-7,
                 'Q4':2.65e-7},
             2: {'S0':3.68e-7,
                 'S1':4.98e-7,
                 'S2':5.60e-7,
                 'S3':5.77e-7,
                 'S4':5.57e-7,},
             3: {'S0':3.88e-7,
                 'S1':5.14e-7,
                 'S2':5.63e-7,
                 'S3':5.63e-7,
                 'S4':5.22e-7,
                 'S5':4.50e-7,}
             }

reverse_transdiff = {v:k for k,v in transdiff.items()}

# def aval(v,ju,jl):
#     """
#     Lookup table for Einstein-A value as a function of
#     vibrational level, upper/lower J level
#     Values from: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
#     """
#     try:
#         return aval_dict[v][reverse_transdiff[ju-jl]+str(jl)]
#     except KeyError:
#         return None


def get_h2_tbl():
    from astroquery.hitran import Hitran
    tbl = Hitran.query_lines(molecule_number=45,
                               isotopologue_number=1,
                               min_frequency=(30*u.um).to(u.THz, u.spectral()),
                               max_frequency=(0.9*u.um).to(u.THz, u.spectral()))
    tbl['branch'] = [x[5] for x in tbl['local_lower_quanta']]
    tbl['jl'] = [int(x.split()[1].strip('q')) for x in tbl['local_lower_quanta']]
    tbl['ju'] = [LL + transdiff[BB] for LL, BB in zip(tbl['jl'], tbl['branch'])]
    tbl['wl'] = (tbl['nu']*u.cm**-1).to(u.um, u.spectral())
    tbl['el_K'] = (tbl['elower']*u.cm**-1 * constants.c * constants.h / constants.k_B).to(u.K)
    tbl['global_lower_quanta'] = tbl['global_lower_quanta'].astype('int')
    tbl['global_upper_quanta'] = tbl['global_upper_quanta'].astype('int')
    return tbl


def aval(vu, vl, ju, jl):
    tbl = get_h2_tbl()
    return tbl['a'][(tbl['jl'] == jl) & (tbl['ju'] == ju) &
            (tbl['global_lower_quanta'].astype('int') == vl) &
            (tbl['global_upper_quanta'].astype('int') == vu)][0] * u.s**-1


def pylevel_population(temperature, orthopararatio, jmax, vmax):

    para = 1/(1.+orthopararatio)
    ortho = 1 - para

    # Compute thermal level populations
    lp = {(V,J):
          (2*J+1)*np.exp(-h2level_energy(V,J)/(kb*temperature))*
          (para if J % 2 == 0 else ortho)
          for J in range(0,jmax)
          for V in range(0,vmax)}
    total_population = sum(list(lp.values()))
    lp = {k:v/total_population for k,v in lp.items()}
    return lp


def emission_per_atom(temperature, orthopararatio, jmax=14, vmax=6):
    """
    Compute the emission per H2 atom in the rovibrational states assuming
    thermal level populations given a temperature and an ortho/para ratio of H2
    """
    lp = level_population(temperature, orthoparatio, jmax, vmax)

    nphotons = {(V,J): aval(V,J,J-2) * lp[(V,J)]
                for J in range(2,14)
                for V in range(1,6)
                if aval(V,J,J-2) is not None}

    return nphotons

def emission_per_atom_oneline(temperature, orthopararatio, vu, ju, jl, jmax=14, vmax=6):
    """
    Compute the emission per H2 atom in the rovibrational states assuming
    thermal level populations given a temperature and an ortho/para ratio of H2
    """

    lp = level_population(temperature, orthopararatio, jmax, vmax)

    nphotons = aval(vu,ju,jl) * lp[(vu,ju)]

    return nphotons

def vectorized_emission_per_atom_oneline(temperatures, orthopararatio, vu, ju,
                                         jl, jmax=14, vmax=6):
    return np.array([emission_per_atom_oneline(t, orthopararatio, vu, ju, jl,
                                               jmax=jmax, vmax=vmax)
                     for t in temperatures])


