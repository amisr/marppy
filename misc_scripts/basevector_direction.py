# basevector_direction.py

import numpy as np
from marppy import Marp
from apexpy import Apex, ApexHeightError
import pymap3d as pm

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs


def basevector_direction():

    glat, glon = np.meshgrid(np.arange(-90., 90., 20.), np.arange(-180., 180., 30.))
    glat = glat.flatten()
    glon = glon.flatten()
    galt = np.full(glat.shape, 300.)

    marp = Marp(date=2019, lam0=80., phi0=45., geo=True)
    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat, glon, galt)

    x, y, z = pm.geodetic2ecef(glat, glon, galt*1000.)
    e1x, e1y, e1z = pm.enu2uvw(e1[0], e1[1], e1[2], glat, glon)
    e2x, e2y, e2z = pm.enu2uvw(e2[0], e2[1], e2[2], glat, glon)
    e3x, e3y, e3z = pm.enu2uvw(e3[0], e3[1], e3[2], glat, glon)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(x, y, z)
    ax.quiver(x, y, z, e1x*1000000., e1y*1000000., e1z*1000000., color='b')
    ax.quiver(x, y, z, e2x*1000000., e2y*1000000., e2z*1000000., color='g')
    ax.quiver(x, y, z, e3x*1000000., e3y*1000000., e3z*1000000., color='r')

    glat, glon, _ = marp.marp2geo(0., 0., 300.)
    x, y, z = pm.geodetic2ecef(glat, glon, galt)
    ax.scatter(x, y, z, marker='^', s=100.)

    plt.show()

def meridians_apex(basevectors=False):

    apexlat = np.arange(-75., 76., 30.)
    apexlon = np.arange(0., 360., 60.)
    galt = 300.

    apex = Apex(date=2019)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')

    # plot magnetic latitude lines
    alon = np.arange(0., 360., 2.)
    for alat in apexlat:
        glat, glon, _ = apex.apex2geo(alat, alon, galt)

        x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
        # x, y, z = pm.geodetic2ecef(alat, alon, np.full(alat.shape, galt*1000.))
        ax.plot(x, y, z)

    # # plot magnetic longitude lines
    # alat = np.arange(-90., 91., 1.)
    # for alon in apexlon:
    #     glat, glon, _ = apex.apex2geo(alat, alon, galt)
    #
    #     x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
    #     # x, y, z = pm.geodetic2ecef(alat, alon, np.full(alat.shape, galt*1000.))
    #     ax.plot(x, y, z)

    if basevectors:
        alat, alon = np.meshgrid(apexlat, apexlon)
        glat, glon, _ = apex.apex2geo(alat.flatten(), alon.flatten(), galt)
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = apex.basevectors_apex(glat, glon, galt)

        # glat = alat.flatten()
        # glon = alon.flatten()
        # d1 = np.broadcast_to([1,0,0],(len(glat),3)).T
        # d2 = np.broadcast_to([0,1,0],(len(glat),3)).T
        # d3 = np.broadcast_to([0,0,1],(len(glat),3)).T
        # e1 = np.broadcast_to([1,0,0],(len(glat),3)).T
        # e2 = np.broadcast_to([0,1,0],(len(glat),3)).T
        # e3 = np.broadcast_to([0,0,1],(len(glat),3)).T


        x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
        e1x, e1y, e1z = pm.enu2uvw(e1[0], e1[1], e1[2], glat, glon)
        e2x, e2y, e2z = pm.enu2uvw(e2[0], e2[1], e2[2], glat, glon)
        e3x, e3y, e3z = pm.enu2uvw(e3[0], e3[1], e3[2], glat, glon)
        d1x, d1y, d1z = pm.enu2uvw(d1[0], d1[1], d1[2], glat, glon)
        d2x, d2y, d2z = pm.enu2uvw(d2[0], d2[1], d2[2], glat, glon)
        d3x, d3y, d3z = pm.enu2uvw(d3[0], d3[1], d3[2], glat, glon)

        ax.scatter(x, y, z)
        # ax.quiver(x, y, z, e1x*1000000., e1y*1000000., e1z*1000000., color='b')
        # ax.quiver(x, y, z, e2x*1000000., e2y*1000000., e2z*1000000., color='g')
        # ax.quiver(x, y, z, e3x*1000000., e3y*1000000., e3z*1000000., color='r')
        ax.quiver(x, y, z, d1x*1000000., d1y*1000000., d1z*1000000., color='b')
        ax.quiver(x, y, z, d2x*1000000., d2y*1000000., d2z*1000000., color='g')
        ax.quiver(x, y, z, d3x*1000000., d3y*1000000., d3z*1000000., color='r')

    plt.show()

def meridians_marp(basevectors=False):

    marplat = np.arange(-80., 81., 20.)
    marplon = np.arange(0., 360., 30.)
    galt = 0.

    # marp = Marp(date=2019, lam0=83.7, phi0=-26.1)
    marp = Marp(date=2019, lam0=30., phi0=0.)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')

    # plot magnetic latitude lines
    mlon = np.arange(0., 360., 2.)
    for mlat in marplat:
        glat, glon, _ = marp.marp2geo(mlat, mlon, galt)
        # glat, glon = marp.marp2apex(mlat, mlon)

        x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
        ax.plot(x, y, z, linewidth=0.5)

    # # plot magnetic longitude lines
    # mlat = np.arange(-90., 91., 1.)
    # for mlon in marplon:
    #     glat, glon, _ = marp.marp2geo(mlat, mlon, galt)
    #     # glat, glon = marp.marp2apex(mlat, mlon)
    #
    #     x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
    #     ax.plot(x, y, z, linewidth=0.5)

    if basevectors:
        mlat, mlon = np.meshgrid(marplat, marplon)
        glat, glon, _ = marp.marp2geo(mlat.flatten(), mlon.flatten(), galt)
        d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat, glon, galt)
        # glat, glon = marp.marp2apex(mlat.flatten(), mlon.flatten())
        # d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat, glon, galt, coords='apex')

        x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
        e1x, e1y, e1z = pm.enu2uvw(e1[0], e1[1], e1[2], glat, glon)
        e2x, e2y, e2z = pm.enu2uvw(e2[0], e2[1], e2[2], glat, glon)
        e3x, e3y, e3z = pm.enu2uvw(e3[0], e3[1], e3[2], glat, glon)
        d1x, d1y, d1z = pm.enu2uvw(d1[0], d1[1], d1[2], glat, glon)
        d2x, d2y, d2z = pm.enu2uvw(d2[0], d2[1], d2[2], glat, glon)
        d3x, d3y, d3z = pm.enu2uvw(d3[0], d3[1], d3[2], glat, glon)

        ax.scatter(x, y, z)
        # ax.quiver(x, y, z, e1x*1000000., e1y*1000000., e1z*1000000., color='b')
        # ax.quiver(x, y, z, e2x*1000000., e2y*1000000., e2z*1000000., color='g')
        # ax.quiver(x, y, z, e3x*1000000., e3y*1000000., e3z*1000000., color='r')
        ax.quiver(x, y, z, d1x*1000000., d1y*1000000., d1z*1000000., color='b')
        ax.quiver(x, y, z, d2x*1000000., d2y*1000000., d2z*1000000., color='g')
        ax.quiver(x, y, z, d3x*1000000., d3y*1000000., d3z*1000000., color='r')

    glat, glon, _ = marp.marp2geo(0., 0., galt)
    x, y, z = pm.geodetic2ecef(glat, glon, galt)
    ax.scatter(x, y, z, marker='^', s=100.)

    x, y, z = pm.geodetic2ecef(90., 0., galt)
    ax.scatter(x, y, z, marker='o', s=100.)

    plt.show()

def zoom():
    mlat, mlon = np.meshgrid(np.arange(-5., 6., 1.), np.arange(-10., 10., 2.))
    mlat = mlat.flatten()
    mlon = mlon.flatten()
    galt = 0.

    marp = Marp(date=2019, lam0=75., phi0=-91., alt=300, geo=True)
    # marp = Marp(date=2019, lam0=0., phi0=0.)
    glat, glon, _ = marp.marp2geo(mlat, mlon, galt)
    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat, glon, galt)
    # glat, glon, _ = marp.apex2geo(mlat, mlon, galt)
    # _,_,_,_,_,_, d1, d2, d3, e1, e2, e3 = marp.basevectors_apex(glat, glon, galt)

    x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
    e1x, e1y, e1z = pm.enu2uvw(e1[0], e1[1], e1[2], glat, glon)
    e2x, e2y, e2z = pm.enu2uvw(e2[0], e2[1], e2[2], glat, glon)
    e3x, e3y, e3z = pm.enu2uvw(e3[0], e3[1], e3[2], glat, glon)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(x, y, z)
    ax.quiver(x, y, z, e1x*100000., e1y*100000., e1z*100000., color='b')
    ax.quiver(x, y, z, e2x*100000., e2y*100000., e2z*100000., color='g')
    ax.quiver(x, y, z, e3x*100000., e3y*100000., e3z*100000., color='r')

    plt.show()

def horiz_vert_components():

    glatitude = np.arange(-90., 91., 5.)
    glongitude = np.arange(-180., 181., 10.)
    glat, glon = np.meshgrid(glatitude, glongitude)

    # apex = Apex(date=2019)
    # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = apex.basevectors_apex(glat.flatten(), glon.flatten(), 300.)
    marp = Marp(date=2019, lam0=40., phi0=0.)
    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 300.)

    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(3,3)

    for i, (name, e) in enumerate([('e1',e1), ('e2',e2), ('e3',e3)]):

        v = e[2,:].reshape(glat.shape)
        h = np.sqrt(e[0,:]**2+e[1,:]**2).reshape(glat.shape)

        ax = plt.subplot(gs[0,i], projection=ccrs.Mollweide())
        c = ax.contourf(glongitude, glatitude, h.T, vmin=-1., vmax=1., cmap=plt.get_cmap('seismic'), transform=ccrs.PlateCarree())
        ax.set_title('{} Horizontal'.format(name))
        plt.colorbar(c)

        ax = plt.subplot(gs[1,i], projection=ccrs.Mollweide())
        c = ax.contourf(glongitude, glatitude, v.T, vmin=-1., vmax=1., cmap=plt.get_cmap('seismic'), transform=ccrs.PlateCarree())
        ax.set_title('{} Vertical'.format(name))
        plt.colorbar(c)

        ax = plt.subplot(gs[2,i], projection=ccrs.Mollweide())
        c = ax.contourf(glongitude, glatitude, np.abs(v.T)/np.abs(h.T), transform=ccrs.PlateCarree())
        ax.set_title('{} Horizontal/Vertical'.format(name))
        plt.colorbar(c)

    # plt.savefig('reference.png')
    plt.show()

def covariant_direction():

    risr_glat = 74.7
    risr_glon = -94.9
    # risr_glat = 0.0
    # risr_glon = 0.0

    # proj = ccrs.Orthographic(central_longitude=risr_glon, central_latitude=risr_glat)
    proj = ccrs.RotatedPole(pole_longitude=risr_glon, pole_latitude=-10., central_rotated_longitude=risr_glon)

    # glat, glon = np.meshgrid(np.arange(70., 90., 5.), np.arange(0., 360., 20.))
    # galt = 0.
    # e1 = np.array([np.full((glat.shape),1.),np.full((glat.shape),0.),np.full((glat.shape),0.)])
    # e2 = np.array([np.full((glat.shape),0.),np.full((glat.shape),1.),np.full((glat.shape),0.)])
    # e3 = np.array([np.full((glat.shape),0.),np.full((glat.shape),0.),np.full((glat.shape),1.)])

    # alat, alon = np.meshgrid(np.arange(70., 90., 5.), np.arange(0., 360., 20.))
    galt = 0.
    alon = np.arange(0., 360., 20.)
    alat = np.full(alon.shape, 80.)
    apex = Apex(date=2019)
    glat, glon, _ = apex.apex2geo(alat, alon, galt)
    _, _, _, _, _, _, d1, d2, d3, e1, e2, e3 = apex.basevectors_apex(glat, glon, galt)
    # _, _, _, _, _, _, d1, d2, d3, e1, e2, e3 = apex.basevectors_apex(glat.ravel(), glon.ravel(), galt)
    # d1 = d1.reshape((3,alat.shape[0],alat.shape[1]))
    # d2 = d2.reshape((3,alat.shape[0],alat.shape[1]))
    # d3 = d3.reshape((3,alat.shape[0],alat.shape[1]))
    # e1 = e1.reshape((3,alat.shape[0],alat.shape[1]))
    # e2 = e2.reshape((3,alat.shape[0],alat.shape[1]))
    # e3 = e3.reshape((3,alat.shape[0],alat.shape[1]))

    # marp = Marp(date=2019, lam0=75., phi0=-91., alt=300, geo=True)
    # lam0 = 83.7
    # phi0 = -26.1
    # lam0 = 0.0
    # phi0 = 0.0
    # marp = Marp(date=2019, lam0=lam0, phi0=phi0, coords='apex')
    # apex = Apex(date=2019)
    # # glat, glon, _ = marp.marp2geo(mlat.flatten(), mlon.flatten(), galt)
    # # d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat, glon, galt)
    # # glat, glon, _ = apex.apex2geo(mlat.flatten(), mlon.flatten(), galt)
    # # _, _, _, _, _, _, d1, d2, d3, e1, e2, e3 = apex.basevectors_apex(glat, glon, galt)
    # mlat, mlon = apex.geo2apex(glat, glon, galt)
    # e1 = e1.reshape((3,mlat.shape[0],mlat.shape[1]))
    # e2 = e2.reshape((3,mlat.shape[0],mlat.shape[1]))
    # glat = glat.reshape(mlat.shape)
    # glon = glon.reshape(mlat.shape)
    print(glat.shape)

    # e1x, e1y, e1z = pm.enu2uvw(e1[0], e1[1], e1[2], glat, glon)
    # e2x, e2y, e2z = pm.enu2uvw(e2[0], e2[1], e2[2], glat, glon)
    # e3x, e3y, e3z = pm.enu2uvw(e3[0], e3[1], e3[2], glat, glon)
    #
    # x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, galt*1000.))
    #
    # e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z, x, y, z, glat, glon = [v.reshape(mlat.shape) for v in [e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z, x, y, z, glat, glon]]
    #
    #
    # dpe = np.diff(mlon,axis=0,append=np.nan)
    # dpn = np.diff(mlon,axis=1,append=np.nan)
    # dpz = np.diff(z,axis=0,append=np.nan)
    # print(diffx.shape, diffy.shape, diffz.shape)
    # dle = np.diff(mlat,axis=0,append=np.nan)
    # dln = np.diff(mlat,axis=1,append=np.nan)
    # dlz = np.diff(z,axis=1,append=np.nan)
    #
    # e1dotdp = (e1x*dpx + e1y*dpy + e1z*dpz)/(np.sqrt(e1x**2+e1y**2+e1z**2)*np.sqrt(dpx**2+dpy**2+dpz**2))
    # e2dotdl = (e2x*dlx + e2y*dly + e2z*dlz)/(np.sqrt(e2x**2+e2y**2+e2z**2)*np.sqrt(dlx**2+dly**2+dlz**2))
    # # print(x.shape, dot.shape)

    fig = plt.figure(figsize=(15,5))
    gs = gridspec.GridSpec(1,2)

    # ax = fig.add_subplot(111,projection='3d')

    # ax.scatter(x, y, z, c=dot.flatten())
    # ax.quiver(x, y, z, e1x*100000., e1y*100000., e1z*100000., color='b')
    # ax.quiver(x, y, z, e2x*100000., e2y*100000., e2z*100000., color='g')
    # ax.quiver(x, y, z, e3x*100000., e3y*100000., e3z*100000., color='r')
    # ax.quiver(x, y, z, diffx, diffy, diffz)

    # # ax = plt.subplot(gs[0],projection=ccrs.Mollweide())
    # # ax = plt.subplot(gs[0],projection=ccrs.LambertConformal(central_longitude=-94.,central_latitude=74.))
    # ax = plt.subplot(gs[0],projection=ccrs.Orthographic(central_longitude=-94.,central_latitude=74.))
    # ax.coastlines()
    # c = ax.scatter(glon, glat, c=e1dotdp, s=1, transform=ccrs.Geodetic())
    # # c = ax.contourf(glon, glat, e1dotdp, 50, transform=ccrs.PlateCarree())
    # ax.scatter(-94., 74., s=10, marker='^', color='red', transform=ccrs.Geodetic())
    # ax.set_title(r'$\hat{{e_1}}\cdot\Delta\phi$, $\lambda_0$ = {}, $\phi_0$ = {}'.format(lam0, phi0))
    # plt.colorbar(c)

    # ax = plt.subplot(gs[0],projection=ccrs.Mollweide())
    # ax = plt.subplot(gs[0],projection=ccrs.LambertConformal(central_longitude=-94.,central_latitude=74.))
    # ax = plt.subplot(gs[0],projection=ccrs.Orthographic(central_longitude=-94.,central_latitude=74.))
    ax = plt.subplot(gs[0],projection=proj)
    ax.coastlines()

    # Geodetic
    c = ax.scatter(glon, glat, s=20, transform=ccrs.Geodetic())
    ax.quiver(glon, glat, e1[0], e1[1], color='red', transform=ccrs.PlateCarree())
    # ax.quiver(glon, glat, e2[0], e2[1], color='blue', transform=ccrs.PlateCarree())

    # ax.quiver(glon, glat, e1[0], e1[1], color='red', transform=ccrs.PlateCarree())
    # ax.quiver(glon, glat, e2[0], e2[1], color='blue', transform=ccrs.PlateCarree())
    # ax.quiver(glon, glat, dpe, dpn, color='red', transform=ccrs.PlateCarree())
    # ax.quiver(glon, glat, e2[0], e2[1], color='blue', transform=ccrs.PlateCarree())
    # c = ax.contourf(glon, glat, e1dotdp, 50, transform=ccrs.PlateCarree())
    ax.scatter(risr_glon, risr_glat, s=50, marker='^', color='orange', transform=ccrs.Geodetic())
    # ax.set_title(r'$\hat{{e_1}}\cdot\Delta\phi$, $\lambda_0$ = {}, $\phi_0$ = {}'.format(lam0, phi0))
    plt.colorbar(c)

    # # ax = plt.subplot(gs[1],projection=ccrs.Mollweide())
    # # ax = plt.subplot(gs[1],projection=ccrs.LambertConformal(central_longitude=-94.,central_latitude=74.))
    # ax = plt.subplot(gs[1],projection=ccrs.Orthographic(central_longitude=-94.,central_latitude=74.))
    # ax.coastlines()
    # c = ax.scatter(glon, glat, c=e2dotdl, s=1, transform=ccrs.Geodetic())
    # # c = ax.contourf(glon, glat, e2dotdl, 50, transform=ccrs.PlateCarree())
    # ax.scatter(-94., 74., s=10, marker='^', color='red', transform=ccrs.Geodetic())
    # ax.set_title(r'$\hat{{e_2}}\cdot\Delta\lambda$, $\lambda_0$ = {}, $\phi_0$ = {}'.format(lam0, phi0))
    # plt.colorbar(c)

    plt.show()

def check_scaling():
    glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(0., 361., 2.))
    # mlat = mlat.flatten()
    # mlon = mlon.flatten()
    galt = 300.

    # marp = Marp(date=2019, lam0=75., phi0=-91., alt=300, geo=True)
    lam0 = 83.7
    phi0 = -26.1
    marp = Marp(date=2019, lam0=10., phi0=0.)
    # glat, glon, _ = marp.marp2geo(mlat.flatten(), mlon.flatten(), galt)
    d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), galt)
    # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = marp.basevectors_apex(glat.flatten(), glon.flatten(), galt)
    # c = np.einsum('i...,i...->...',np.cross(d1.T,d2.T).T,d3).reshape(glat.shape)
    c = np.linalg.norm(np.cross(d1.T,d2.T).T)/np.linalg.norm(e3)

    # np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',np.cross(d1.T,d2.T).T,d3), np.ones(glat.flatten().shape), decimal=2)


    print(c[np.abs(c-1.0)>0.0001])
    print(c.shape, d1.shape, glat.shape)

    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(1,1)

    # for i, (name, e) in enumerate([('e1',e1), ('e2',e2), ('e3',e3)]):

        # v = e[2,:].reshape(glat.shape)
        # h = np.sqrt(e[0,:]**2+e[1,:]**2).reshape(glat.shape)

    ax = plt.subplot(gs[0], projection=ccrs.Mollweide())
    c = ax.scatter(glon, glat, c=c, s=1., cmap=plt.get_cmap('jet'), transform=ccrs.PlateCarree())
    # ax.set_title('{} Horizontal'.format(name))
    plt.colorbar(c)

    plt.show()




def main():
    # basevector_direction()
    # meridians_apex(basevectors=True)
    # meridians_marp(basevectors=True)
    # zoom()
    horiz_vert_components()
    # covariant_direction()
    # check_scaling()

if __name__=='__main__':
    main()
