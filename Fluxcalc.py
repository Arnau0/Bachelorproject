import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

c = 29979245800  # cm/s
h = 6.62507015 / (10**27)  # cm^2 g s^-1
k_B = 1.380649 / (10**16)  # cm^2 g s^-2 K^-1
R_sun = 6.96 * (10**10)  # cm
pc = 3.08567758 * (10**18)  # cm
au = 1.49597871 * (10**13)  # cm

# example values from Kavanagh et al. (2020)
exT = 10**6  # K
exnu = 10**6  # Hz = 100 MHz
exn_0 = 10**8  # cm^-3


def abs_coeff_r(T, nu, n_0, r):
    """calculates the absorption coefficient alpha_nu that is dependent on r through the relation n = n_0*r^-2.
    For fully ionized hydrogen, n_e = n_i & Z=1
    T = temperature
    nu = frequency
    n_0 = initial number density profile
    r = distance to absorption (from observer?) => variable
    """
    C = 3.692 * (10**8)  # the  number in front of the expression
    BB = 1 - np.exp(
        -(h * nu) / (k_B * T)
    )  # term that looks like a term in the black body spectrum
    g = 10.6 + 1.9 * np.log10(T) + 1.26 * np.log10(nu)  # Gaunt factor
    return (C * BB * g * (n_0**2)) / (np.sqrt(T) * (nu**3) * 4 * (r**4 / R_sun**4))


L_grid = 500
lmax = 100 * R_sun
x_axis = np.linspace(-lmax, lmax, L_grid)
y_axis = np.linspace(0, lmax, int(L_grid / 2))

# note the index for wind_coordinate_system: x = [:,0] & y = [:,1]
x, y = np.meshgrid(x_axis, y_axis, indexing="ij")
r = (x**2 + y**2) ** 0.5  # 2d grid with r vector
mask = (r < R_sun) | ((x > 0) & (abs(y) < R_sun))


def I_tau(tau_max, nu):  # the intensity given a maximum optical depth
    B = ((2 * h * nu**3) / (c**2)) * (np.exp((h * nu) / (k_B * exT)) - 1) ** (-1)
    return B * (1 - np.exp(-tau_max))


frequency = np.linspace(10**6, (10**10), L_grid)

dx = x_axis[1] - x_axis[0]
dy = y_axis[1] - y_axis[0]
dA = 2 * np.pi * y_axis * dy

Flux_nu = []
for nu in frequency:

    # absorption coefficient
    alpha_nu = abs_coeff_r(exT, nu, exn_0, r)
    alpha_nu[mask] = 0

    tau_max = np.sum(alpha_nu, axis=0) * dx

    # intensity
    I_nu = I_tau(tau_max, nu)

    # flux
    F = np.sum(I_nu * dA) / (pc**2)
    Flux_nu.append(F)

plt.plot(frequency, Flux_nu)
plt.xlabel("frequency")
plt.ylabel("Flux")
plt.xscale("log")
plt.yscale("log")
plt.savefig("Flux-frequency.png")
