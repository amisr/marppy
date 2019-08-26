# coord_map.py

import numpy as np
from apexpy import Apex
from marp import Marp

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

def Richmond1995_fig3():
    A = Apex(1995)

    lat, lon = np.meshgrid(np.arange(-90.,90.,5.),np.arange(-180.,180.,5.))

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(lat.ravel(), lon.ravel(), 110., coords='geo')

    d1xd2 = np.cross(d1.T, d2.T).T
    D = np.linalg.norm(d1xd2, axis=0).reshape(lat.shape)

    d12 = np.einsum('i...,i...->...', d1, d1).reshape(lat.shape)/D
    d1d2 = np.einsum('i...,i...->...', d1, d2).reshape(lat.shape)/D
    d22 = np.einsum('i...,i...->...', d2, d2).reshape(lat.shape)/D

    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2,2)

    for i, val in enumerate([D,d12,d1d2,d22]):
        ax = plt.subplot(gs[i], projection=ccrs.PlateCarree())
        ax.coastlines()
        c = ax.contourf(lon, lat, val, 30, cmap=plt.get_cmap('jet'), transform=ccrs.PlateCarree())

        plt.colorbar(c)
    plt.show()


def coord_map():

    lat, lon = np.meshgrid(np.arange(-90.,91.,5.),np.arange(0.,361.,5.))
    lat = lat.flatten()
    lon = lon.flatten()


    A = Apex(2019)
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(lat, lon, 0., coords='apex')
    d1xd2 = np.cross(d1.T, d2.T).T
    D = np.linalg.norm(d1xd2, axis=0)
    d12 = np.einsum('i...,i...->...', d1, d1)/D
    d1d2 = np.einsum('i...,i...->...', d1, d2)/D
    d22 = np.einsum('i...,i...->...', d2, d2)/D
    apex = [lat, D, d12, d1d2, d22]
    labels = ['MLAT', r'$D$', r'$d_1^2/D$', r'$d_1\cdot d_2/D$', r'$d_2^2/D$']


    marp = []
    marp_lat = []
    marp_lon = []
    marp_title = []
    for (lat0, lon0) in [(0.,60.), (30.,60.), (60.,60.), (90.,60.)]:
        M = Marp(date=2019, lam0=lat0, phi0=lon0)

        rplat, rplon = M.apex2marp(lat, lon)
        marp_d1, marp_d2, marp_d3, marp_e1, marp_e2, marp_e3 = M.basevectors_marp(lat, lon, 0., coords='apex')
        marp_d1xd2 = np.cross(marp_d1.T, marp_d2.T).T
        # marp_D = np.linalg.norm(marp_d1xd2, axis=0)
        marp_D = np.linalg.norm(marp_e3, axis=0)
        marp_d12 = np.einsum('i...,i...->...', marp_d1, marp_d1)/marp_D
        marp_d1d2 = np.einsum('i...,i...->...', marp_d1, marp_d2)/marp_D
        marp_d22 = np.einsum('i...,i...->...', marp_d2, marp_d2)/marp_D
        marp.append([lat, marp_D, marp_d12, marp_d1d2, marp_d22])
        marp_lat.append(rplat)
        marp_lon.append(rplon)
        marp_title.append('({}, {})'.format(lat0,lon0))


    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(len(apex), len(marp)+1)
    gs.update(left=0.05,right=0.92,top=0.99,bottom=0.01,hspace=0.01,wspace=0.1)

    for i in range(len(apex)):
        ax = plt.subplot(gs[i,0], projection=ccrs.Mollweide())
        c = ax.scatter(lon, lat, c=apex[i], cmap=plt.get_cmap('jet'), s=2.5, transform=ccrs.PlateCarree())
        ax.set_title('Apex')
        ax.text(-0.25, 0.5, labels[i], transform=ax.transAxes)

        for j in range(len(marp)):
            ax = plt.subplot(gs[i,j+1], projection=ccrs.Mollweide())
            ax.scatter(marp_lon[j], marp_lat[j], c=marp[j][i], vmin=min(apex[i]), vmax=max(apex[i]), cmap=plt.get_cmap('jet'), s=2.5, transform=ccrs.PlateCarree())
            ax.set_title(marp_title[j])

        axpos = ax.get_position() # get the original position 
        cax = fig.add_axes([axpos.x1 + 0.01, axpos.y0,  0.02, axpos.y1-axpos.y0])
        plt.colorbar(c, cax=cax)

    plt.savefig('coord_map.png')
    plt.show()

def main():
    # Richmond1995_fig3()
    coord_map()

if __name__=='__main__':
    main()