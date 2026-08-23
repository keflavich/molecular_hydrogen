#!/usr/bin/env python
"""
Regenerate ``molecular_hydrogen/data/h2_level_energies.ecsv``.

HITRAN tabulates, for every line, the lower state energy and the transition
wavenumber, along with the (v, J) assignment of both levels.  Every X-state
level that appears in the line list can therefore be recovered as either
``elower`` (when it is the lower level of a transition) or ``elower + nu``
(when it is the upper level).  Most levels are reached by many transitions, so
averaging the determinations also provides a consistency check.

This is run manually; the output is committed so that the package works
offline.
"""
import numpy as np
from astropy import units as u
from astropy.table import Table
from astroquery.hitran import Hitran

TRANSDIFF = {'S': 2, 'R': 1, 'Q': 0, 'P': -1, 'O': -2}


def main(output='molecular_hydrogen/data/h2_level_energies.ecsv'):
    tbl = Hitran.query_lines(molecule_number=45, isotopologue_number=1,
                             min_frequency=(1e6*u.um).to(u.THz, u.spectral()),
                             max_frequency=(0.1*u.um).to(u.THz, u.spectral()))

    vl = np.array([int(x) for x in tbl['global_lower_quanta']])
    vu = np.array([int(x) for x in tbl['global_upper_quanta']])
    jl = np.array([int(x.split()[1].rstrip('qm'))
                   for x in tbl['local_lower_quanta']])
    branch = np.array([x[5] for x in tbl['local_lower_quanta']])
    ju = np.array([j + TRANSDIFF[b] for j, b in zip(jl, branch)])
    elower = np.array(tbl['elower'], dtype='float64')
    nu = np.array(tbl['nu'], dtype='float64')

    determinations = {}
    for v, j, e in zip(vl, jl, elower):
        determinations.setdefault((v, j), []).append(e)
    for v, j, e, n in zip(vu, ju, elower, nu):
        determinations.setdefault((v, j), []).append(e + n)

    rows = []
    for (v, j), values in sorted(determinations.items()):
        # HITRAN stores the lower state energy to 4 decimal places; the
        # float32 columns astroquery returns carry noise beyond that
        rows.append((v, j, round(float(np.mean(values)), 4),
                     round(float(np.ptp(values)), 4), len(values)))

    out = Table(rows=rows, names=['V', 'J', 'energy', 'spread', 'ndet'])
    out['energy'].unit = u.cm**-1
    out['energy'].description = 'term value above the (V=0, J=0) level'
    out['spread'].unit = u.cm**-1
    out['spread'].description = 'range of the independent determinations'
    out['ndet'].description = 'number of transitions the level appears in'
    out.meta['source'] = ('HITRAN2020 (molecule 45, isotopologue 1), '
                          'retrieved via astroquery.hitran')
    out.meta['generated_by'] = 'tools/make_level_energies.py'
    out.write(output, overwrite=True)

    print("%i levels, V=0-%i, J=0-%i" % (len(out), out['V'].max(), out['J'].max()))
    print("worst disagreement between determinations: %.4g cm^-1"
          % out['spread'].max())


if __name__ == '__main__':
    main()
