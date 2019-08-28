# test.py

import unittest
import numpy as np
from marp import Marp

class TestMarp(unittest.TestCase):

    def test_init(self):
        marp = Marp(date=2019, lam0=60., phi0=30.)
        self.assertIsInstance(marp, Marp)

    def test_apex2marp_norotation(self):
        # validate that if no rotation is set (lam0=0 & phi0=0), marp coordinates equal apex coordinates
        marp = Marp(date=2019, lam0=0., phi0=0.)
        mlat, mlon = marp.apex2marp(60., 30.)

        self.assertAlmostEqual(mlat, 60.)
        self.assertAlmostEqual(mlon, 30.)

    def test_apex2marp_center(self):
        # validate that the center point becomes (0,0) in marp
        marp = Marp(date=2019, lam0=60., phi0=30.)
        mlat, mlon = marp.apex2marp(60., 30.)

        self.assertAlmostEqual(mlat, 0.)
        self.assertAlmostEqual(mlon, 0.)

    def test_marp2apex_norotation(self):
        # validate that if no rotation is set (lam0=0 & phi0=0), apex coordinates equal marp coordinates
        marp = Marp(date=2019, lam0=0., phi0=0.)
        alat, alon = marp.marp2apex(60., 30.)

        self.assertAlmostEqual(alat, 60.)
        self.assertAlmostEqual(alon, 30.)

    def test_marp2apex_center(self):
        # validate that marp (0,0) becomes the center point in apesx
        marp = Marp(date=2019, lam0=60., phi0=30.)
        alat, alon = marp.marp2apex(0., 0.)

        self.assertAlmostEqual(alat, 60.)
        self.assertAlmostEqual(alon, 30.)

    def test_geo2marp(self):
        # validate that tranforming from geodetic to marp and back to geodetic returns the original coordinates
        marp = Marp(date=2019, lam0=60., phi0=30.)
        mlat, mlon = marp.geo2marp(50., 10., 300.)
        glat, glon, _ = marp.marp2geo(mlat, mlon, 300.)

        self.assertAlmostEqual(glat, 50., places=4)
        self.assertAlmostEqual(glon, 10., places=4)

    def test_marp2geo(self):
        # validate that transforming from marp to geodetic and back to marp returns the original coordinates
        marp = Marp(date=2019, lam0=60., phi0=30.)
        glat, glon, _ = marp.marp2geo(50., 10., 300.)
        mlat, mlon = marp.geo2marp(glat, glon, 300.)

        self.assertAlmostEqual(mlat, 50., places=4)
        self.assertAlmostEqual(mlon, 10., places=4)

    def test_basevectors_marp_reciprocal(self):
        # validate that the basevectors are reciprocal, i.e. di*ej=delta{ij} (Richmond 1995, eqn 3.18)
        marp = Marp(date=2019, lam0=60., phi0=30.)
        d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(50., 10., 300.)

        self.assertAlmostEqual(np.dot(d1,e1), 1., places=5)
        self.assertAlmostEqual(np.dot(d1,e2), 0., places=5)
        self.assertAlmostEqual(np.dot(d1,e3), 0., places=5)
        self.assertAlmostEqual(np.dot(d2,e1), 0., places=5)
        self.assertAlmostEqual(np.dot(d2,e2), 1., places=5)
        self.assertAlmostEqual(np.dot(d2,e3), 0., places=5)
        self.assertAlmostEqual(np.dot(d3,e1), 0., places=5)
        self.assertAlmostEqual(np.dot(d3,e2), 0., places=5)
        self.assertAlmostEqual(np.dot(d3,e3), 1., places=5)

    def test_basevectors_marp_cross(self):
        # validate the property d1xd2*d3 = e1xe2*e3 = 1 (Richmond 1995, eqn 3.17)
        marp = Marp(date=2019, lam0=60., phi0=30.)
        d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(50., 10., 300.)

        self.assertAlmostEqual(np.dot(np.cross(d1, d2), d3), 1., places=5)
        self.assertAlmostEqual(np.dot(np.cross(e1, e2), e3), 1., places=5)

    def test_basevectors_marp_de_relation(self):
        # validate that the relationships betwee d and e shown in Richmond 1995, eqn 3.16 
        #   and Laundal and Richmond 2016, eqn 57-59 hold true
        marp = Marp(date=2019, lam0=60., phi0=30.)
        d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(50., 10., 300.)

        np.testing.assert_almost_equal(np.cross(e2, e3), d1)
        np.testing.assert_almost_equal(np.cross(e3, e1), d2)
        np.testing.assert_almost_equal(np.cross(e1, e2), d3)

        np.testing.assert_almost_equal(np.cross(d2, d3), e1)
        np.testing.assert_almost_equal(np.cross(d3, d1), e2)
        np.testing.assert_almost_equal(np.cross(d1, d2), e3)

    def test_mapping(self):
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

        self.assertAlmostEqual(Ed1A, Ed1B, places=2)
        self.assertAlmostEqual(Ed2A, Ed2B, places=2)
        self.assertAlmostEqual(Ve1A, Ve1B, places=2)
        self.assertAlmostEqual(Ve2A, Ve2B, places=2)


if __name__=='__main__':
    unittest.main()