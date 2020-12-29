# derivative_check.py
# This is a short script inteded to validate the partial derivaties used to calculate MARP base vectors.
# The script uses sympy to symbolically calculate the relevant partial derivaties and compare them with
#   the expressions copied from marp.py.
# If all expressions are correct, all fields of the results printed to screen will be zero.

from sympy import *

def contravariant_derivatives():

    x, y, z, zr, yr, zr, l, p, lr, pr, l0, p0 = symbols('x y z xr yr zr l p lr pr l0 p0')

    xr = cos(l0)*cos(l)*cos(p-p0) + sin(l0)*sin(l)
    yr = cos(l)*sin(p-p0)
    zr = -sin(l0)*cos(l)*cos(p-p0) + cos(l0)*sin(l)

    pr = atan2(yr,xr)
    lr = asin(zr)

    # partial derivatives hard-coded in marp.py
    dprdp = cos(l)*(cos(l0)*cos(l)+sin(l0)*sin(l)*cos(p-p0))/(xr**2+yr**2)
    dprdl = -sin(l0)*sin(p-p0)/(xr**2+yr**2)
    dlrdp = sin(l0)*cos(l)*sin(p-p0)/sqrt(1-zr**2)
    dlrdl = (sin(l0)*sin(l)*cos(p-p0)+cos(l0)*cos(l))/sqrt(1-zr**2)

    print('CONTRAVARIANT DERIVATIVES:')
    print('dprdp:', simplify(diff(pr, p)-dprdp))
    print('dprdl:', simplify(diff(pr, l)-dprdl))
    print('dlrdp:', simplify(diff(lr, p)-dlrdp))
    print('dlrdl:', simplify(diff(lr, l)-dlrdl))

def covariant_derivatives():

    x, y, z, zr, yr, zr, l, p, lr, pr, l0, p0 = symbols('x y z xr yr zr l p lr pr l0 p0')

    x = cos(p0)*cos(l0)*cos(lr)*cos(pr) - cos(p0)*sin(l0)*sin(lr) - sin(p0)*cos(lr)*sin(pr)
    y = sin(p0)*cos(l0)*cos(lr)*cos(pr) - sin(p0)*sin(l0)*sin(lr) + cos(p0)*cos(lr)*sin(pr)
    z = sin(l0)*cos(lr)*cos(pr) + cos(l0)*sin(lr)

    p = atan2(y,x)
    l = asin(z)

    # partial derivatives hard-coded in marp.py
    dpdpr = cos(lr)*(cos(l0)*cos(lr)-sin(l0)*sin(lr)*cos(pr))/(x**2+y**2)
    dpdlr = sin(l0)*sin(pr)/(x**2+y**2)
    dldpr = -sin(l0)*cos(lr)*sin(pr)/sqrt(1-z**2)
    dldlr = (-sin(l0)*sin(lr)*cos(pr)+cos(l0)*cos(lr))/sqrt(1-z**2)

    print('COVARIANT DERIVATIVES')
    print('dpdpr:', simplify(diff(p, pr)-dpdpr))
    print('dpdlr:', simplify(diff(p, lr)-dpdlr))
    print('dldpr:', simplify(diff(l, pr)-dldpr))
    print('dldlr:', simplify(diff(l, lr)-dldlr))

def main():
    contravariant_derivatives()
    covariant_derivatives()

if __name__=='__main__':
    main()
