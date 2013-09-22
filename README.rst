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

   $ pip install -e git@github.com:keflavich/molecular_hydrogen.git

Usage
-----

.. code-block:: python

   >>> from molecular_hydrogen import h2
   >>> h2.linename_to_restwl('1-0 Q(1)')
   2.406594067745623
   >>> import astropy.units as u
   >>> h2.linename_to_restwl('S(0) 0-0')*u.um
   <Quantity 28.2206857627 um>
   

.. _Download ZIP: https://github.com/keflavich/molecular_hydrogen/zipball/master
.. _Download TAR: https://github.com/keflavich/molecular_hydrogen/tarball/master
.. _View on Github: https://github.com/keflavich/molecular_hydrogen/
.. _Documentation: https://github.com/keflavich/molecular_hydrogen/README.rst


.. image:: https://d2weczhvl823v0.cloudfront.net/keflavich/molecular_hydrogen/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

