# test_marppy.py

import pytest
import numpy as np
from marppy import Marp
from apexpy import Apex


def test_apex2marp_norotation():
    # validate that if no rotation is set (lam0=0 & phi0=0), marp coordinates equal apex coordinates
    marp = Marp(date=2019, lam0=0., phi0=0.)
    alat, alon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-178., 182., 2.))

    mlat, mlon = marp.apex2marp(alat, alon)

    np.testing.assert_array_almost_equal(mlat, alat)
    np.testing.assert_array_almost_equal(mlon, alon)


def test_marp2apex_norotation():
    # validate that if no rotation is set (lam0=0 & phi0=0), apex coordinates equal marp coordinates
    marp = Marp(date=2019, lam0=0., phi0=0.)
    mlat, mlon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-178., 182., 2.))

    alat, alon = marp.apex2marp(mlat, mlon)

    np.testing.assert_array_almost_equal(alat, mlat)
    np.testing.assert_array_almost_equal(alon, mlon)


def test_apex2marp_center():
    # validate that the center point becomes (0,0) in marp
    marp = Marp(date=2019, lam0=60., phi0=30.)
    mlat, mlon = marp.apex2marp(60., 30.)

    np.testing.assert_almost_equal(mlat, 0.)
    np.testing.assert_almost_equal(mlon, 0.)


def test_marp2apex_center():
    # validate that marp (0,0) becomes the center point in apex
    marp = Marp(date=2019, lam0=60., phi0=30.)
    alat, alon = marp.marp2apex(0., 0.)

    np.testing.assert_almost_equal(alat, 60.)
    np.testing.assert_almost_equal(alon, 30.)


def test_apex2marp():
    # validate that trasforming from apex to marp and back to apex returns the original coordinates
    marp = Marp(date=2019, lam0=60., phi0=30.)

    alat0, alon0 = np.meshgrid(np.arange(-89., 90., 1.), np.arange(-178., 180., 2.))

    mlat, mlon = marp.apex2marp(alat0, alon0)
    alat, alon = marp.marp2apex(mlat, mlon)

    np.testing.assert_array_almost_equal(alat, alat0, decimal=3)
    np.testing.assert_array_almost_equal(alon, alon0, decimal=3)


def test_marp2apex():
    # validate that trasforming from marp to apex and back to marp returns the original coordinates
    marp = Marp(date=2019, lam0=60., phi0=30.)

    mlat0, mlon0 = np.meshgrid(np.arange(-89., 90., 1.), np.arange(-178., 180., 2.))

    alat, alon = marp.marp2apex(mlat0, mlon0)
    mlat, mlon = marp.apex2marp(alat, alon)

    np.testing.assert_array_almost_equal(mlat, mlat0, decimal=3)
    np.testing.assert_array_almost_equal(mlon, mlon0, decimal=3)


def test_geo2marp():
    # validate that tranforming from geodetic to marp and back to geodetic returns the original coordinates
    marp = Marp(date=2019, lam0=60., phi0=30.)

    glat0, glon0 = np.meshgrid(np.arange(-89., 90., 1.), np.arange(-178., 180., 2.))

    mlat, mlon = marp.geo2marp(glat0, glon0, 0.)
    glat, glon, _ = marp.marp2geo(mlat, mlon, 0.)

    np.testing.assert_array_almost_equal(glat, glat0, decimal=3)
    np.testing.assert_array_almost_equal(glon, glon0, decimal=3)


def test_marp2geo():
    # validate that transforming from marp to geodetic and back to marp returns the original coordinates
    marp = Marp(date=2019, lam0=60., phi0=30.)

    mlat0, mlon0 = np.meshgrid(np.arange(-89., 90., 1.), np.arange(-178., 180., 2.))

    glat, glon, _ = marp.marp2geo(mlat0, mlon0, 0.)
    mlat, mlon = marp.geo2marp(glat, glon, 0.)

    np.testing.assert_array_almost_equal(mlat, mlat0, decimal=3)
    np.testing.assert_array_almost_equal(mlon, mlon0, decimal=3)


def test_basevectors_norotation():
    # validate that if no rotation is set (lam0=0 & phi0=0), marp base vectors equal apex base vectors
    apex = Apex(date=2019)
    marp = Marp(date=2019, lam0=0., phi0=0.)

    # glat, glon = np.meshgrid(np.arange(-89., 90., 2.), np.arange(-179., 180., 4.))
    glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-180., 181., 2.))

    _, _, _, _, _, _, Ad1, Ad2, Ad3, Ae1, Ae2, Ae3 = apex.basevectors_apex(glat.flatten(), glon.flatten(), 0.)
    Md1, Md2, Md3, Me1, Me2, Me3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)

    np.testing.assert_array_almost_equal(Ad1, Md1)
    np.testing.assert_array_almost_equal(Ad2, Md2)
    np.testing.assert_array_almost_equal(Ad3, Md3)
    np.testing.assert_array_almost_equal(Ae1, Me1)
    np.testing.assert_array_almost_equal(Ae2, Me2)
    np.testing.assert_array_almost_equal(Ae3, Me3)


def test_basevectors_reciprocal():
    # validate that the basevectors are reciprocal, i.e. di*ej=delta{ij} (Richmond 1995, eqn 3.18)
    marp = Marp(date=2019, lam0=60., phi0=30.)
    glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-180., 181., 2.))

    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)

    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d1,e1), np.ones(glat.flatten().shape), decimal=3)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d1,e2), np.zeros(glat.flatten().shape), decimal=2)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d1,e3), np.zeros(glat.flatten().shape), decimal=2)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d2,e1), np.zeros(glat.flatten().shape), decimal=2)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d2,e2), np.ones(glat.flatten().shape), decimal=3)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d2,e3), np.zeros(glat.flatten().shape), decimal=2)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d3,e1), np.zeros(glat.flatten().shape), decimal=2)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d3,e2), np.zeros(glat.flatten().shape), decimal=2)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d3,e3), np.ones(glat.flatten().shape), decimal=3)


# def test_general_co_contra_relationship(self):
#     # validate the general relationship ei = ejxek/ejxek*ei (Arfken - Exercise 4.3.1)
#     # This should hold even if base vectors are not properly scaled
#     marp = Marp(date=2019, lam0=60., phi0=30.)
#     glat, glon = np.meshgrid(np.arange(-89., 90., 2.), np.arange(-179., 180., 4.))
#
#     d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)
#
#     eps3 = np.cross(e1.T,e2.T).T/np.einsum('i...,i...->...',np.cross(e1.T,e2.T).T,e3)
#     np.testing.assert_array_almost_equal(eps3, d3, decimal=4)
#     eps1 = np.cross(e2.T,e3.T).T/np.einsum('i...,i...->...',np.cross(e2.T,e3.T).T,e1)
#     np.testing.assert_array_almost_equal(eps1, d1, decimal=4)
#     eps2 = np.cross(e3.T,e1.T).T/np.einsum('i...,i...->...',np.cross(e3.T,e1.T).T,e2)
#     np.testing.assert_array_almost_equal(eps2, d2, decimal=4)


def test_basevectors_scaled():
    # validate that the base vectors are properly scaled such that d1xd2*d3 = e1xe2*e3 = 1 (Richmond 1995, eqn 3.17)
    # this is a special property of Apex base vectors
    marp = Marp(date=2019, lam0=60., phi0=30.)
    # glat, glon = np.meshgrid(np.arange(50., 90., 2.), np.arange(-179., 180., 4.))
    glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-180., 181., 2.))

    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)

    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',np.cross(d1.T,d2.T).T,d3), np.ones(glat.flatten().shape), decimal=4)
    np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',np.cross(e1.T,e2.T).T,e3), np.ones(glat.flatten().shape), decimal=4)


def test_basevectors_parallel2B():
    # validate that d3 and e3 are parallel
    marp = Marp(date=2019, lam0=60., phi0=30.)
    glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-180., 181., 2.))

    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)

    d3u = d3/np.linalg.norm(d3,axis=0)
    e3u = e3/np.linalg.norm(e3,axis=0)

    np.testing.assert_array_almost_equal(d3u, e3u, decimal=4)


def test_basevectors_construction():
    # validate that the relationships between d and e shown in Richmond 1995, eqn 3.16
    #   and Laundal and Richmond 2016, eqn 57-59 hold true
    marp = Marp(date=2019, lam0=60., phi0=30.)
    glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-180., 181., 2.))

    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)

    np.testing.assert_almost_equal(np.cross(e2.T, e3.T).T, d1, decimal=2)
    np.testing.assert_almost_equal(np.cross(e3.T, e1.T).T, d2, decimal=2)
    np.testing.assert_almost_equal(np.cross(e1.T, e2.T).T, d3, decimal=2)

    np.testing.assert_almost_equal(np.cross(d2.T, d3.T).T, e1, decimal=2)
    np.testing.assert_almost_equal(np.cross(d3.T, d1.T).T, e2, decimal=2)
    np.testing.assert_almost_equal(np.cross(d1.T, d2.T).T, e3, decimal=2)


def test_mapping():
    # validate that for an Electric field mapped to different altitudes, Ed1 and Ed2 are constant and
    #   for a Velocity field maped to different altitudes, Ve1 and Ve2 are constant
    marp = Marp(date=2019, lam0=60., phi0=30.)

    alat = 70.
    alon = 60.
    altA = 100.
    altB = 4000.

    # set an arbitrary electric field and use standard apexpy functions to map it
    EA = np.array([300.,100.,0.])
    EB = marp.map_E_to_height(alat, alon, altA, altB, EA)
    # set an arbitrary velocity field and use standard apexpy functions to map it
    VA = np.array([200.,400.,0.])
    VB = marp.map_V_to_height(alat, alon, altA, altB, VA)

    # calculate basevectors and components at altitude A
    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(alat, alon, altA, coords='apex')
    Ed1A = np.dot(EA, e1)
    Ed2A = np.dot(EA, e2)
    Ve1A = np.dot(VA, d1)
    Ve2A = np.dot(VA, d2)

    # calculate basevectors and components at altitude B
    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(alat, alon, altB, coords='apex')
    Ed1B = np.dot(EB, e1)
    Ed2B = np.dot(EB, e2)
    Ve1B = np.dot(VB, d1)
    Ve2B = np.dot(VB, d2)

    np.testing.assert_almost_equal(Ed1A, Ed1B, decimal=3)
    np.testing.assert_almost_equal(Ed2A, Ed2B, decimal=3)
    np.testing.assert_almost_equal(Ve1A, Ve1B, decimal=3)
    np.testing.assert_almost_equal(Ve2A, Ve2B, decimal=3)
