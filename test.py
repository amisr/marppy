# test.py

import numpy as np
from marp import Marp

def test_conversion():

    rp = Marp(date=2019, lam0=60., phi0=30.)

    lat = np.arange(0., 90., 1.)
    lon = np.arange(0., 180., 2.)
    alt = np.full(lat.shape, 0.)
    print(lat, lon, alt)

    alat, alon = rp.geo2apex(lat, lon, alt)
    mlat, mlon = rp.apex2marp(alat, alon)
    alat, alon = rp.marp2apex(mlat, mlon)
    glat, glon, _ = rp.apex2geo(alat, alon, alt)
    print(glat-lat, glon-lon)

    mlat, mlon = rp.geo2marp(lat, lon, alt)
    glat, glon, _ = rp.marp2geo(mlat, mlon, alt)
    print(glat-lat, glon-lon)

def test_basevectors():

    rp = Marp(date=2019, lam0=60., phi0=30.)

    lat = np.arange(0., 90., 1.)
    lon = np.arange(0., 180., 2.)
    alt = np.full(lat.shape, 0.)

    d1, d2, d3, e1, e2, e3 = rp.basevectors_marp(lat, lon, alt, coords='geo')

    d1xd2 = np.cross(d1.T, d2.T).T

    # print('d3-d1xd2', d3-d1xd2/np.linalg.norm(d1xd2, axis=0)**2)
    print('e3-d1xd2', e3-d1xd2)

    print('d1.e1', np.mean(np.einsum('i...,i...->...', d1, e1)))
    print('d1.e2', np.mean(np.einsum('i...,i...->...', d1, e2)))
    print('d1.e3', np.mean(np.einsum('i...,i...->...', d1, e3)))
    print('d2.e1', np.mean(np.einsum('i...,i...->...', d2, e1)))
    print('d2.e2', np.mean(np.einsum('i...,i...->...', d2, e2)))
    print('d2.e3', np.mean(np.einsum('i...,i...->...', d2, e3)))
    print('d3.e1', np.mean(np.einsum('i...,i...->...', d3, e1)))
    print('d3.e2', np.mean(np.einsum('i...,i...->...', d3, e2)))
    print('d3.e3', np.mean(np.einsum('i...,i...->...', d3, e3)))

def test_mapping():

    rp = Marp(date=2019, lam0=60., phi0=30.)

    alat = 70.
    alon = 60.
    EfieldA = np.array([500.,500.,0.])
    EfieldB = rp.map_E_to_height(alat, alon, 100., 4000., EfieldA)
    print(EfieldA, EfieldB)

    # evaluate at MAGNETIC lat/lons
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = rp.basevectors_apex(alat, alon, 100., coords='apex')
    Ed1A = np.einsum('i,i->', EfieldA, e1)
    Ed2A = np.einsum('i,i->', EfieldA, e2)

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = rp.basevectors_apex(alat, alon, 4000., coords='apex')
    Ed1B = np.einsum('i,i->', EfieldB, e1)
    Ed2B = np.einsum('i,i->', EfieldB, e2)

    print('APEX:', Ed1A, Ed2A, Ed1B, Ed2B)


    d1, d2, d3, e1, e2, e3 = rp.basevectors_marp(alat, alon, 100., coords='apex')
    Ed1A = np.einsum('i,i->', EfieldA, e1)
    Ed2A = np.einsum('i,i->', EfieldA, e2)

    d1, d2, d3, e1, e2, e3 = rp.basevectors_marp(alat, alon, 4000., coords='apex')
    Ed1B = np.einsum('i,i->', EfieldB, e1)
    Ed2B = np.einsum('i,i->', EfieldB, e2)

    print('MARP:', Ed1A, Ed2A, Ed1B, Ed2B)


def main():
    # test_conversion()
    test_basevectors()
    test_mapping()

if __name__=='__main__':
    main()