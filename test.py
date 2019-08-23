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

def main():
    test_conversion()

if __name__=='__main__':
    main()