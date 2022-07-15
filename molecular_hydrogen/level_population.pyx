from astropy import constants
import pyximport
pyximport.install()
from h2level_energy import h2level_energy
import numpy as np

cdef float h,c,kb,e,We,Be,WeXe,De,Ae

h = constants.h.cgs.value #6.626068e-27
c = constants.c.cgs.value #2.99792e10
kb = constants.k_B.cgs.value
e = constants.e.esu.value #4.803e-12

We=4401.21
Be=60.853
WeXe=121.33
De=0.0471
Ae=3.062


def level_population(float temperature, float orthopararatio, int jmax,
                     int vmax):

    cdef float para,ortho,total_population,v
    cdef int V,J

    para = 1./(1.+orthopararatio)
    ortho = 1. - para

    # Compute thermal level populations
    cdef float kbt = (kb*temperature)
    if kbt == 0.0:
        # can't have div by zero...
        return {(V,J):np.nan
                for J in range(0,jmax)
                for V in range(0,vmax)}

    else:
        lp = {(V,J):
                    (2*J+1)*np.exp(-h2level_energy(V,J)/kbt)*
                    (para if J % 2 == 0 else ortho)
                    for J in range(0,jmax)
                    for V in range(0,vmax)}
        vals = np.array(list(lp.values()), dtype='float')
        total_population = vals.sum()
        if total_population == 0:
            return lp
        else:
            lp = {k:v/total_population for k,v in lp.items()}
            return lp
