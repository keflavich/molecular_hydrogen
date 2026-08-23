`Documentation`_ |  `View on Github`_ |  `Download ZIP`_  |  `Download TAR`_  

Molecular Hydrogen Properties
=============================

A python toolkit to compute molecular hydrogen properties: rest wavelength,
Einstein A coefficients, etc.  

Right now, supports only the rovibrational transitions and is strongly targeted
at the near-IR regime.

Why?
----
The NIST database:
http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740
http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1000
has none of these values, just the tools to produce them.

There are other pages:
http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
that have nice lists of values, but they aren't easily machine readable.


Finally, the original source of these values is very difficult to find.

Installation
------------

.. code-block:: bash

   $ pip install -e git+git@github.com:keflavich/molecular_hydrogen.git#egg=molecular_hydrogen

``numpy`` and ``astropy`` are required.  ``astroquery`` is needed for the
HITRAN lookups (Einstein A coefficients, full line lists) and ``cython`` is
optional -- it speeds up the level population calculation about 8x, and the
package falls back to a pure-python implementation when it is unavailable.

Usage
-----

Rest wavelengths are returned as ``astropy`` quantities:

.. code-block:: python

   >>> from molecular_hydrogen import h2
   >>> h2.linename_to_restwl('1-0 Q(1)')
   <Quantity 2.40659407 micron>
   >>> h2.linename_to_restwl('S(0) 0-0')
   <Quantity 28.22068576 micron>
   >>> h2.restwl(vu=1, vl=0, ju=3, jl=1)
   <Quantity 2.12183131 micron>

Einstein A coefficients come from HITRAN:

.. code-block:: python

   >>> h2.aval(vu=1, vl=0, ju=3, jl=1)              # doctest: +SKIP
   <Quantity 3.467e-07 1 / s>

The full HITRAN line list, with the quantum numbers decoded, is available as a
table.  For example, every transition between 1 and 6 microns:

.. code-block:: python

   >>> tbl = h2.get_h2_tbl(min_wavelength=1*u.um,
   ...                     max_wavelength=6*u.um)    # doctest: +SKIP
   >>> tbl['linename', 'wl', 'a', 'el_K'][:3]        # doctest: +SKIP

Note that HITRAN lists the Q-branch lines twice, as electric quadrupole
(``transition_type == 'q'``) and magnetic dipole (``'m'``) transitions.  The
quadrupole lines are the astrophysically interesting ones and are what ``aval``
returns by default.

Data sources
------------

Rest wavelengths are computed from the measured X-state term values of
Dabrowski (1984), Can. J. Phys. 62, 1639, covering v = 0-5 and J = 0-13.
Einstein A coefficients and the extended line list come from HITRAN
(molecule 45).

Tests
-----

.. code-block:: bash

   $ pytest                 # offline tests
   $ pytest --remote-data   # also cross-check against HITRAN


.. _Download ZIP: https://github.com/keflavich/molecular_hydrogen/zipball/master
.. _Download TAR: https://github.com/keflavich/molecular_hydrogen/tarball/master
.. _View on Github: https://github.com/keflavich/molecular_hydrogen/
.. _Documentation: https://github.com/keflavich/molecular_hydrogen/README.rst


.. image:: https://d2weczhvl823v0.cloudfront.net/keflavich/molecular_hydrogen/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

