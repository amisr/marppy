# reciprocal_check.py

from marppy import Marp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

marp = Marp(date=2019, lam0=80., phi0=30.)
glat, glon = np.meshgrid(np.arange(-90., 91., 1.), np.arange(-180., 181., 2.))

d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat.flatten(), glon.flatten(), 0.)
# d1, d2, d3, e1, e2, e3 = marp.basevectors_marp(glat, glon, 0.)

fig = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(3,3)

for i, d in enumerate([d1, d2, d3]):
    for j, e in enumerate([e1, e2, e3]):
        # if i==j:
        #     vmin = 1.-(1.e-7)
        #     vmax = 1.+(1.e-7)
        # else:
        vmin = 1-1.e-7
        vmax = 1+1.e-7
        r = np.einsum('i...,i...->...',d,e).reshape(glat.shape)
        ax = fig.add_subplot(gs[i,j])
        if i==j:
            c = ax.pcolormesh(glon, glat, r, vmin=1-1.e-6, vmax=1+1.e-6, cmap='coolwarm')
        else:
            c = ax.pcolormesh(glon, glat, r, vmin=-1.e-6, vmax=1.e-6, cmap='coolwarm')
        fig.colorbar(c)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d1,e1), np.ones(glat.flatten().shape), decimal=3)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d1,e2), np.zeros(glat.flatten().shape), decimal=2)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d1,e3), np.zeros(glat.flatten().shape), decimal=2)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d2,e1), np.zeros(glat.flatten().shape), decimal=2)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d2,e2), np.ones(glat.flatten().shape), decimal=3)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d2,e3), np.zeros(glat.flatten().shape), decimal=2)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d3,e1), np.zeros(glat.flatten().shape), decimal=2)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d3,e2), np.zeros(glat.flatten().shape), decimal=2)
# np.testing.assert_array_almost_equal(np.einsum('i...,i...->...',d3,e3), np.ones(glat.flatten().shape), decimal=3)

plt.show()
