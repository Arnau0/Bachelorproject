import numpy as np
import math as math
import matplotlib.pyplot as plt

npts = 1000
r = np.logspace(0, 2, npts)
z = 3 * np.sin(r) + r**0.3

coord = np.linspace(-50, 50, 42)
x_grid, y_grid = np.meshgrid(coord, coord, indexing="ij")
r_grid = (x_grid**2 + y_grid**2) ** 0.5

# method 1: numpy interp
z_grid1 = np.interp(r_grid, r, z)


# method 2: scipy interp1d
# "This class is considered legacy and will no longer receive updates.
# This could also mean it will be removed in future SciPy versions."
from scipy.interpolate import interp1d

z_interp = interp1d(r, z)
z_grid2 = z_interp(r_grid)


from matplotlib.colors import LogNorm

fig = plt.figure()

fig.add_subplot(2, 1, 1)
plt.imshow(z_grid1, norm=LogNorm())
plt.colorbar()

fig.add_subplot(2, 1, 2)
plt.imshow(z_grid2, norm=LogNorm())
plt.colorbar()
plt.show()
