from astropy import constants

cdef float h,c,k,kb,e,We,Be,WeXe,De,Ae

h = constants.h.cgs.value #6.626068e-27
c = constants.c.cgs.value #2.99792e10
kb = constants.k_B.cgs.value
e = constants.e.esu.value #4.803e-12

We=4401.21
Be=60.853
WeXe=121.33
De=0.0471
Ae=3.062

def h2level_energy(int V,int J):
    """ Returns the theoretical level energy as a function of the
    vibrational (V) and rotational (J) state of the molecule.
    
    Constants are from NIST:
    http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1000#Diatomic
    (see the bottom of the table)

    Returns a value in ergs
    """

    return (h * c *
            (We*(V+0.5) + 
             Be*(J*(J+1.)) - 
             WeXe*(V+.5)**2. - 
             De*J**2.*(J+1)**2. - 
             Ae*(V+.5)*(J+1)*J))

