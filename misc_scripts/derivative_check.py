# derivative_check.py
# This is a short script inteded to validate the partial derivaties used to calculate MARP base vectors.
# The script uses sympy to symbolically calculate the relevant partial derivaties and compare them with
#   the expressions copied from marp.py.
# If all expressions are correct, all fields of the results printed to screen will be zero.

import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def contravariant_derivatives():

    # x, y, z, zr, yr, zr, l, p, lr, pr, l0, p0 = symbols('x y z xr yr zr l p lr pr l0 p0')
    #
    # xr = cos(l0)*cos(l)*cos(p-p0) + sin(l0)*sin(l)
    # yr = cos(l)*sin(p-p0)
    # zr = -sin(l0)*cos(l)*cos(p-p0) + cos(l0)*sin(l)
    #
    # pr = atan2(yr,xr)
    # lr = asin(zr)
    #
    # # partial derivatives hard-coded in marp.py
    # dprdp = cos(l)*(cos(l0)*cos(l)+sin(l0)*sin(l)*cos(p-p0))/(xr**2+yr**2)
    # dprdl = -sin(l0)*sin(p-p0)/(xr**2+yr**2)
    # dlrdp = sin(l0)*cos(l)*sin(p-p0)/sqrt(1-zr**2)
    # dlrdl = (sin(l0)*sin(l)*cos(p-p0)+cos(l0)*cos(l))/sqrt(1-zr**2)
    #
    # print('CONTRAVARIANT DERIVATIVES:')
    # print('dprdp:', simplify(diff(pr, p)-dprdp))
    # print('dprdl:', simplify(diff(pr, l)-dprdl))
    # print('dlrdp:', simplify(diff(lr, p)-dlrdp))
    # print('dlrdl:', simplify(diff(lr, l)-dlrdl))


    # Plot derivitative caluculated analytically and numerically to check correctness of analytic solution

    lam0 = 60.*np.pi/180.
    phi0 = 30.*np.pi/180.
    tau0 = 20.*np.pi/180.

    pA, lA = np.meshgrid(np.arange(0., 360., 0.5)*np.pi/180., np.arange(-90., 90., 0.5)*np.pi/180.)

    xA = np.cos(lA)*np.cos(pA)
    yA = np.cos(lA)*np.sin(pA)
    zA = np.sin(lA)
    rA = np.array([xA, yA, zA]).T

    Rtau = np.array([[np.cos(tau0), np.sin(tau0), 0.], [-np.sin(tau0), np.cos(tau0), 0.], [0., 0., 1.]])
    Rlam = np.array([[np.sin(lam0), 0., -np.cos(lam0)], [0., 1., 0.], [np.cos(lam0), 0., np.sin(lam0)]])
    Rphi = np.array([[np.cos(phi0), np.sin(phi0), 0.], [-np.sin(phi0), np.cos(phi0), 0.], [0., 0., 1.]])
    # print(rr)
    R = np.einsum('ij,jk,kl->il', Rtau, Rlam, Rphi)

    rM = np.einsum('ij,...j->...i', R, rA).T

    xM = rM[0]
    yM = rM[1]
    zM = rM[2]

    pM = np.arctan2(yM, xM)
    lM = np.arcsin(zM)


    P1 = np.array([[-np.sin(pA)*np.cos(lA), -np.cos(pA)*np.sin(lA)], [np.cos(pA)*np.cos(lA), -np.sin(pA)*np.sin(lA)], [np.zeros(lA.shape), np.cos(lA)]])
    P2 = np.array([[-yM/(xM**2+yM**2), xM/(xM**2+yM**2), np.zeros(xM.shape)], [np.zeros(xM.shape), np.zeros(xM.shape), 1./np.sqrt(1-zM**2)]])
    dMdA = np.einsum('ij...,jk,kl...->il...', P2, R, P1)

    print(pA.shape, lA.shape, dMdA.shape)

    dpMdpA = dMdA[0,0]
    dpMdlA = dMdA[0,1]
    dlMdpA = dMdA[1,0]
    dlMdlA = dMdA[1,1]

    num_dpMdlA, num_dpMdpA = np.gradient(pM, lA[:,0], pA[0,:])
    num_dlMdlA, num_dlMdpA = np.gradient(lM, lA[:,0], pA[0,:])
    print(num_dpMdpA.shape)
    # plt.plot(pA, lM)
    # plt.plot(pA, np.gradient(lM, pA))
    # plt.plot(pA, dlMdpA, linestyle=':')
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(4,3)

    ax = fig.add_subplot(gs[0,0])
    ax.pcolormesh(pA, lA, dpMdpA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[0,1])
    ax.pcolormesh(pA, lA, num_dpMdpA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[0,2])
    ax.pcolormesh(pA, lA, dpMdpA-num_dpMdpA, vmin=-1e-2, vmax=1e-2)

    ax = fig.add_subplot(gs[1,0])
    ax.pcolormesh(pA, lA, dpMdlA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[1,1])
    ax.pcolormesh(pA, lA, num_dpMdlA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[1,2])
    ax.pcolormesh(pA, lA, dpMdlA-num_dpMdlA, vmin=-1e-2, vmax=1e-2)

    ax = fig.add_subplot(gs[2,0])
    ax.pcolormesh(pA, lA, dlMdpA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[2,1])
    ax.pcolormesh(pA, lA, num_dlMdpA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[2,2])
    ax.pcolormesh(pA, lA, dlMdpA-num_dlMdpA, vmin=-1e-2, vmax=1e-2)

    ax = fig.add_subplot(gs[3,0])
    ax.pcolormesh(pA, lA, dlMdlA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[3,1])
    ax.pcolormesh(pA, lA, num_dlMdlA, vmin=-np.pi, vmax=np.pi)
    ax = fig.add_subplot(gs[3,2])
    ax.pcolormesh(pA, lA, dlMdlA-num_dlMdlA, vmin=-1e-2, vmax=1e-2)

    plt.show()




def covariant_derivatives():

    x, y, z, zr, yr, zr, l, p, lr, pr, l0, p0 = sym.symbols('x y z xr yr zr l p lr pr l0 p0')

    x = sym.cos(p0)*sym.cos(l0)*sym.cos(lr)*sym.cos(pr) - sym.cos(p0)*sym.sin(l0)*sym.sin(lr) - sym.sin(p0)*sym.cos(lr)*sym.sin(pr)
    y = sym.sin(p0)*sym.cos(l0)*sym.cos(lr)*sym.cos(pr) - sym.sin(p0)*sym.sin(l0)*sym.sin(lr) + sym.cos(p0)*sym.cos(lr)*sym.sin(pr)
    z = sym.sin(l0)*sym.cos(lr)*sym.cos(pr) + sym.cos(l0)*sym.sin(lr)

    p = sym.atan2(y,x)
    l = sym.asin(z)

    # partial derivatives hard-coded in marp.py
    dpdpr = sym.cos(lr)*(sym.cos(l0)*sym.cos(lr)-sym.sin(l0)*sym.sin(lr)*sym.cos(pr))/(x**2+y**2)
    dpdlr = sym.sin(l0)*sym.sin(pr)/(x**2+y**2)
    dldpr = -sym.sin(l0)*sym.cos(lr)*sym.sin(pr)/sym.sqrt(1-z**2)
    dldlr = (sym.cos(l0)*sym.cos(lr)-sym.sin(l0)*sym.sin(lr)*sym.cos(pr))/sym.sqrt(1-z**2)

    print('COVARIANT DERIVATIVES')
    print('dpdpr:', sym.simplify(sym.diff(p, pr)-dpdpr))
    print('dpdlr:', sym.simplify(sym.diff(p, lr)-dpdlr))
    print('dldpr:', sym.simplify(sym.diff(l, pr)-dldpr))
    print('dldlr:', sym.simplify(sym.diff(l, lr)-dldlr))

def main():
    contravariant_derivatives()
    # covariant_derivatives()

if __name__=='__main__':
    main()
