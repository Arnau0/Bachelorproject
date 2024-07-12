# choose M_dot and T
# calculate and plot the flux density

import numpy as np
import math as math
import matplotlib.pyplot as plt
from isothermal_wind import isothermal_wind

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
exT = 2e6  # K
exnu = 1e6  # Hz = 100 MHz
exn_0 = 1e8  # cm^-3
exMdot = 100 * Mdot_sun  # g/s


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


# intensity
def I_tau(tau_max, nu, T):
    B = ((2 * h * nu**3) / (c**2)) * (np.exp((h * nu) / (k_B * T)) - 1) ** (-1)
    return B * (1 - np.exp(-tau_max))


# unit conversion
def toMilliJansky(F):
    return F * int(1e26)


# presets
Mstar = 0.2
Rstar = 0.3
Rmax = 100
L_grid = 500

# range of the system
lmax = Rmax

# the array for the y-axis is only for the top half
x_axis = np.linspace(-lmax, lmax, L_grid)
y_axis = np.linspace(0, lmax, int(L_grid / 2))

x, y = np.meshgrid(x_axis, y_axis, indexing="ij")
r_grid = (x**2 + y**2) ** 0.5  # 2d grid with r vector
mask = (r_grid < Rstar) | ((x > 0) & (abs(y) < Rstar))

dx = (x_axis[1] - x_axis[0]) * R_sun
dy = (y_axis[1] - y_axis[0]) * R_sun
dA = 2 * np.pi * y_axis * dy * R_sun

###########################
# varying temperature
###########################
temperature = np.linspace(1e6, 3e6, L_grid)

F_T = []
for T in temperature:

    # massloss rate and frequency fixed to example values
    r, u = isothermal_wind(Mstar, Rstar, T, Rmax, npts=L_grid)

    # number density profile interpolation
    rho = exMdot / (4 * np.pi * (r * Rstar * R_sun) ** 2 * (u * 1e5))
    n = rho / m_p

    n_grid = np.interp(r_grid, r, n)

    alpha = alpha_n(T, exnu, n_grid)
    alpha[mask] = 0

    tau_max = np.sum(alpha, axis=0) * dx

    I = I_tau(tau_max, exnu, T)

    F = np.sum(I * dA) / (pc**2)
    F = toMilliJansky(F)
    F_T.append(F)


###########################
# varying massloss rate
###########################
massloss = np.linspace(50 * Mdot_sun, 500 * Mdot_sun, L_grid)

F_Mdot = []
for Mdot in massloss:

    # temperature and frequency fixed to example values
    r, u = isothermal_wind(Mstar, Rstar, exT, Rmax, npts=L_grid)

    rho = Mdot / (4 * np.pi * (r * Rstar * R_sun) ** 2 * (u * 1e5))
    n = rho / m_p

    n_grid = np.interp(r_grid, r, n)

    # frequency fixed
    alpha = alpha_n(exT, exnu, n_grid)
    alpha[mask] = 0

    tau_max = np.sum(alpha, axis=0) * dx

    I = I_tau(tau_max, exnu, exT)

    F = np.sum(I * dA) / (pc**2)
    F = toMilliJansky(F)
    F_Mdot.append(F)


plt.plot(temperature, F_T)
plt.xlabel("Temperature (K)")
plt.ylabel("Flux density (mJy)")
# plt.savefig("Flux-temperature.png")
plt.show()

plt.figure()
plt.plot(massloss, F_Mdot)
plt.xlabel("Massloss rate (g/s)")
plt.ylabel("Flux density (mJy)")
# plt.savefig("Flux-massloss.png")
plt.show()
