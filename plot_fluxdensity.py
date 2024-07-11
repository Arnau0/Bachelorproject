# choose M_dot and T
# calculate and plot the flux density

import numpy as np
import math as math
import matplotlib.pyplot as plt
import isothermal_wind

# constants in cgs units
c = 29979245800  # cm/s
h = 6.62507015e-27  # cm^2 g s^-1
k_B = 1.380649e-16  # cm^2 g s^-2 K^-1
pc = 3.08567758e18  # cm
au = 1.49597871e13  # cm
m_p = 1.6726e-24

# solar constants
R_sun = 6.955e10  # cm
Mdot_sun = 1.261e12  # g/s
M_sun = 1.988e33  # g

# example values from Kavanagh et al. (2020)
exT = 1e6  # K
exnu = 1e6  # Hz = 100 MHz
exn_0 = 1e8  # cm^-3


def alpha_n(T, nu, n):
    """calculates the absorption coefficient alpha_nu that is dependent on n through the relation n = 2*n_e.
    T = temperature (K)
    nu = frequency (Hz)
    n = number density profile value (cm^-3)
    """
    C = 3.692 * (10**8)  # the  number in front of the expression
    BB = 1 - np.exp(
        -(h * nu) / (k_B * T)
    )  # term that looks like the term found in the blackbody spectrum
    g = 10.6 + 1.9 * np.log10(T) + 1.26 * np.log10(nu)  # Gaunt factor
    return C * BB * g * ((n**2) / (np.sqrt(T) * (nu**3) * 4))


Mstar = 0.2
Rstar = 0.3
T = 2e6
Rmax = 100
Mdot = 100 * Mdot_sun

L_grid = 500

r, u = isothermal_wind.isothermal_wind(Mstar, Rstar, T, Rmax, npts=L_grid)

# range of the system
lmax = Rmax * R_sun

# the array for the y-axis is only for the top half
x_axis = np.linspace(-lmax, lmax, L_grid)
y_axis = np.linspace(0, lmax, int(L_grid / 2))

x, y = np.meshgrid(x_axis, y_axis, indexing="ij")
r_grid = (x**2 + y**2) ** 0.5  # 2d grid with r vector
mask = (r_grid < R_sun) | ((x > 0) & (abs(y) < R_sun))


# number density profile interpolation
rho = Mdot / (4 * np.pi * (r * Rstar * R_sun) ** 2 * (u * 1e5))
n = rho / m_p

n_grid = np.interp(r_grid, r, n)

alpha = alpha_n(exT, exnu, n_grid)

from matplotlib.colors import LogNorm

plt.imshow(n_grid, norm=LogNorm())
plt.colorbar()
plt.show()
