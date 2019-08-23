# marp.py

import numpy as np
from apexpy import Apex
import pymap3d as pm

RE = 6371.
hR = 0.

class Marp(Apex):
    def __init__(self, date=None, refh=0, datafile=None, lam0=0., phi0=0.):

        self.lam0 = lam0
        self.phi0 = phi0

        super(Marp, self).__init__(date=date, refh=refh, datafile=None)

    def apex2marp(self, lam, phi):

        lam = lam*np.pi/180.
        phi = phi*np.pi/180.
        lam0 = self.lam0*np.pi/180.
        phi0 = self.phi0*np.pi/180.

        xr = np.cos(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) + np.cos(phi0)*np.sin(lam0)*np.sin(lam) + np.sin(phi0)*np.cos(lam)*np.sin(phi)
        yr = -np.sin(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) - np.sin(phi0)*np.sin(lam0)*np.sin(lam) + np.cos(phi0)*np.cos(lam)*np.sin(phi)
        zr = -np.sin(lam0)*np.cos(lam)*np.cos(phi) + np.cos(lam0)*np.sin(lam)

        phir = np.arctan2(yr, xr)
        lamr = np.arcsin(zr)

        return lamr*180./np.pi, phir*180./np.pi

    def marp2apex(self, lamr, phir):

        lamr = lamr*np.pi/180.
        phir = phir*np.pi/180.
        lam0 = self.lam0*np.pi/180.
        phi0 = self.phi0*np.pi/180.

        x = np.cos(lam0)*np.cos(phi0)*np.cos(lamr)*np.cos(phir) - np.cos(lam0)*np.sin(phi0)*np.cos(lamr)*np.sin(phir) - np.sin(lam0)*np.sin(lamr)
        y = np.sin(phi0)*np.cos(lamr)*np.cos(phir) + np.cos(phi0)*np.cos(lamr)*np.sin(phir)
        z = np.sin(lam0)*np.cos(phi0)*np.cos(lamr)*np.cos(phir) - np.sin(lam0)*np.sin(phi0)*np.cos(lamr)*np.sin(phir) + np.cos(lam0)*np.sin(lamr)

        phi = np.arctan2(y, x)
        lam = np.arcsin(z)

        return lam*180./np.pi, phi*180./np.pi

    def geo2marp(self, glat, glon, height):
        alat, alon = self.geo2apex(glat, glon, height)
        mlat, mlon = self.apex2marp(alat, alon)
        return mlat, mlon

    def marp2geo(self, mlat, mlon, height):
        alat, alon = self.marp2apex(mlat, mlon)
        glat, glon, err = self.apex2geo(alat, alon, height)
        return glat, glon, err



def base_vectors(lam, phi, lam0, phi0):
    Apx = Apex(2019)
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = Apx.basevectors_apex(lam, phi, 0., coords='apex')
    glat, glon, _ = Apx.apex2geo(lam, phi, 0.)
    u, v, w = pm.enu2uvw(d1[0], d1[1], d1[2], glat, glon)
    d1 = np.array([u, v, w])
    u, v, w = pm.enu2uvw(d2[0], d2[1], d2[2], glat, glon)
    d2 = np.array([u, v, w])

    lamr, phir = ma2marp(lam, phi, lam0, phi0)

    lam = lam*np.pi/180.
    phi = phi*np.pi/180.
    lam0 = lam0*np.pi/180.
    phi0 = phi0*np.pi/180.
    lamr = lamr*np.pi/180.
    phir = phir*np.pi/180.

    # xr = np.cos(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) + np.cos(phi0)*np.sin(lam0)*np.sin(lam) + np.sin(phi0)*np.cos(lam)*np.sin(phi)
    # yr = -np.sin(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) - np.sin(phi0)*np.sin(lam0)*np.sin(lam) + np.cos(phi0)*np.cos(lam)*np.sin(phi)
    # zr = -np.sin(lam0)*np.cos(lam)*np.cos(phi) + np.cos(lam0)*np.sin(lam)

    # phir = np.arctan2(yr, xr)
    # lamr = np.arcsin(zr)

    d1r = 1/np.cos(lamr)*((np.sin(phir+phi0)*np.cos(lam0)*np.cos(lam)*np.sin(phi)+np.cos(phir+phi0)*np.cos(lam)*np.cos(phi))*d1 + (np.sin(phir+phi0)*np.cos(lam0)*np.sin(lam)*np.cos(phi)-np.sin(phir+phi0)*np.sin(lam0)*np.cos(lam)-np.cos(phir+phi0)*np.sin(lam)*np.sin(phi))*d2)
    d2r = 1/np.cos(lamr)*(np.sin(lam0)*np.cos(lam)*np.sin(phi)*d1 + (np.sin(lam0)*np.sin(lam)*np.cos(phi)+np.cos(lam0)*np.cos(lam))*d2)

    d1xd2 = np.cross(d1r.T, d2r.T).T
    D = np.linalg.norm(d1xd2, axis=0)
    d3r = d1xd2/D**2

    # Currently, this returns ECEF components

    return d1r, d2r, d3r

# def base_vectors2(lam, phi, lam0, phi0):
#     Apx = Apex(2019)
#     f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = Apx.basevectors_apex(lam, phi, 0., coords='apex')
#     glat, glon, _ = Apx.apex2geo(lam, phi, 0.)
#     u, v, w = pm.enu2uvw(d1[0], d1[1], d1[2], glat, glon)
#     d1 = np.array([u, v, w])
#     u, v, w = pm.enu2uvw(d2[0], d2[1], d2[2], glat, glon)
#     d2 = np.array([u, v, w])

#     lam = lam*np.pi/180.
#     phi = phi*np.pi/180.
#     lam0 = lam0*np.pi/180.
#     phi0 = phi0*np.pi/180.

#     xr = np.cos(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) + np.cos(phi0)*np.sin(lam0)*np.sin(lam) + np.sin(phi0)*np.cos(lam)*np.sin(phi)
#     yr = -np.sin(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) - np.sin(phi0)*np.sin(lam0)*np.sin(lam) + np.cos(phi0)*np.cos(lam)*np.sin(phi)
#     zr = -np.sin(lam0)*np.cos(lam)*np.cos(phi) + np.cos(lam0)*np.sin(lam)

#     phir = np.arctan2(yr, xr)
#     lamr = np.arcsin(zr)

#     A = np.sqrt(4-3*np.cos(lam)**2)/2.
#     # d1r = 1./((RE+hR)*np.cos(lamr))*(np.sin(phir+phi0)*np.cos(lam0)*(np.sin(phi)*d1-A*np.cos(phi)*d2)+np.cos(phir+phi0)*(np.cos(phi)*d1+A*np.sin(phi0)*d2)+np.sin(phir+phi0)*np.sin(lam0)*A/np.tan(lam)*d2)
#     # d2r = 1./((RE+hR)*np.cos(lamr))*(np.sin(lam0)*(np.sin(phi)*d1-A*np.cos(phi)*d2)-np.cos(lam0)*A/np.tan(lam)*d2)
#     d1r = (np.sin(phir+phi0)*np.cos(lam0)*(np.sin(phi)*d1-A*np.cos(phi)*d2)+np.cos(phir+phi0)*(np.cos(phi)*d1+A*np.sin(phi0)*d2)+np.sin(phir+phi0)*np.sin(lam0)*A/np.tan(lam)*d2)
#     d2r = -2*np.sin(lamr)/(np.sqrt(4-3*np.cos(lamr)**2)*np.cos(lamr))*(np.sin(lam0)*(np.sin(phi)*d1-A*np.cos(phi)*d2)-np.cos(lam0)*A/np.tan(lam)*d2)

#     return d1r, d2r

# def base_vectors(lam, phi, lam0, phi0):

#     A = Apex(2019)

#     lam = lam*np.pi/180.
#     phi = phi*np.pi/180.
#     lam0 = lam0*np.pi/180.
#     phi0 = phi0*np.pi/180.

#     xr = np.cos(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) + np.cos(phi0)*np.sin(lam0)*np.sin(lam) + np.sin(phi0)*np.cos(lam)*np.sin(phi)
#     yr = -np.sin(phi0)*np.cos(lam0)*np.cos(lam)*np.cos(phi) - np.sin(phi0)*np.sin(lam0)*np.sin(lam) + np.cos(phi0)*np.cos(lam)*np.sin(phi)
#     zr = -np.sin(lam0)*np.cos(lam)*np.cos(phi) + np.cos(lam0)*np.sin(lam)

#     phir = np.arctan2(yr, xr)
#     lamr = np.arcsin(zr)

#     dxrdp = -np.cos(phi0)*np.cos(lam0)*np.cos(lam)*np.sin(phi) + np.sin(phi0)*np.cos(lam)*np.cos(phi)
#     dyrdp = np.sin(phi0)*np.cos(lam0)*np.cos(lam)*np.sin(phi) + np.cos(phi0)*np.cos(lam)*np.cos(phi)
#     dzrdp = np.sin(lam0)*np.cos(lam)*np.sin(phi)

#     dxrdl = -np.cos(phi0)*np.cos(lam0)*np.sin(lam)*np.cos(phi) + np.cos(phi0)*np.sin(lam0)*np.cos(lam) - np.sin(phi0)*np.sin(lam)*np.sin(phi)
#     dyrdl = np.sin(phi0)*np.cos(lam0)*np.sin(lam)*np.cos(phi) - np.sin(phi0)*np.sin(lam0)*np.cos(lam) - np.cos(phi0)*np.sin(lam)*np.sin(phi)
#     dzrdl = np.sin(lam0)*np.sin(lam)*np.cos(phi) + np.cos(lam0)*np.cos(lam)


#     dprdp = (dyrdp*xr-yr*dxrdp)/(xr**2+yr**2)
#     dprdl = (dyrdl*xr-yr*dxrdl)/(xr**2+yr**2)

#     dlrdp = dzrdp/np.sqrt(1-zr**2)
#     dlrdl = dzrdl/np.sqrt(1-zr**2)

#     f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(lam*180./np.pi, phi*180./np.pi, 0., coords='apex')
#     glat, glon, _ = A.apex2geo(lam, phi, 0.)
#     u, v, w = pm.enu2uvw(d1[0], d1[1], d1[2], glat, glon)
#     d1 = np.array([u, v, w])
#     u, v, w = pm.enu2uvw(d2[0], d2[1], d2[2], glat, glon)
#     d2 = np.array([u, v, w])

#     gradp = d1/np.cos(lam)
#     gradl = -np.sqrt(4-3*np.cos(lam)**2)*d2/(2*np.sin(lam))

#     gradpr = dprdp*gradp + dprdl*gradl
#     gradlr = dlrdp*gradp + dlrdl*gradl

#     # d1r = np.cos(lamr)*gradpr

#     # d2r = -2*np.sin(lamr)/np.sqrt(4-3*np.cos(lamr)**2)*gradlr

#     d1r = gradpr
#     d2r = gradlr

#     # d1r = np.cos(lamr)/(xr**2+yr**2)*((np.cos(lam)*np.cos(lam0)+np.cos(phi)*np.sin(lam0)*np.sin(lam))*d1 + np.sin(lam0)*np.sin(phi)*np.sqrt(4-3*np.cos(lam)**2)/(2*np.sin(lam))*d2)

#     # d2r = -2*np.sin(lamr)/(np.sqrt(4-3*np.cos(lamr)**2)*np.sqrt(1-zr**2))*(np.sin(lam0)*np.sin(phi)*d1 - (np.sin(lam0)*np.sin(lam)*np.cos(phi)+np.cos(lam0)*np.cos(lam))*np.sqrt(4-3*np.cos(lam)**2)/(2*np.sin(lam))*d2)
#     # print(d1, d2)

#     return d1r, d2r

