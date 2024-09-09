# Compute the mass-loss rate constraint given a series of lowest detected flux densities
# requires csv file with columns: id, mass, radius, frequency, lowest flux density and distance
# returns the given csv file plus constrainted mass-loss rate


import numpy as np
import math as math
import matplotlib.pyplot as plt
from isothermal_wind import isothermal_wind
import pandas as pd

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
exnu = 1e6  # Hz
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


def getLum(F, d):
    return 4 * np.pi * F * d**2


#####################################
# radiostars data
data = pd.read_csv(
    r"C:\Users\Arnau\Documents\Bachelorproject\Python\Radiostars_for_flux.csv",
    delimiter=",",
)

ID = data["ID"].values.astype(str)
Masses = data["Mass"].values.astype(float)  # in solar masses
Radii = data["Radius"].values.astype(float)  # in solar radii
Frequencies = data["Frequency"].values.astype(float)  # in MHz
Flux = data["Flux_I"].values.astype(float)  # in mJy
Distances = data["Distance"].values.astype(float)  # in pc


def Fluxcalc(Rstar, Mstar, nu, d, Rmax=1500, L_grid=500):
    """Calculate the flux density of a star using the isothermal wind model for the density profile.
    Rstar = radius of the star in solar radii
    Mstar = mass of the star in solar masses
    nu = central frequency in Hz
    d = distance in cm
    Rmax = system size in solar radii
    L_grid = number of points on each side

    Returns list of massloss-rates and corresponding flux densities.
    """
    # 500 grid points is quite large if you have ~200 stars to process (about 5 minutes total)

    #####################################
    # system grid

    # the array for the y-axis is only for the top half
    x_axis = np.linspace(-Rmax, Rmax, L_grid)
    y_axis = np.linspace(0, Rmax, int(L_grid / 2))

    x, y = np.meshgrid(x_axis, y_axis, indexing="ij")
    r_grid = (x**2 + y**2) ** 0.5  # 2d grid with r vector
    mask = (r_grid < Rstar) | ((x > 0) & (abs(y) < Rstar))

    dx = (x_axis[1] - x_axis[0]) * R_sun * Rstar
    dy = (y_axis[1] - y_axis[0]) * R_sun * Rstar
    dA = 2 * np.pi * y_axis * dy * R_sun * Rstar
    dA[0] = np.pi * dy**2

    #####################################
    # calculating the flux density

    massloss = np.logspace(0, 7, 100) * Mdot_sun
    # massloss = np.linspace(1 * Mdot_sun, 400000 * Mdot_sun, L_grid)

    # temperature fixed to value: (example value exT=2e6 K)
    r, u = isothermal_wind(Mstar, Rstar, exT, Rmax, npts=L_grid)

    F_Mdot = []
    for Mdot in massloss:

        rho = Mdot / (4 * np.pi * (r * Rstar * R_sun) ** 2 * (u * 1e5))
        n = rho / m_p

        n_grid = np.interp(r_grid, r, n)

        # frequency fixed
        alpha = alpha_n(exT, nu, n_grid)
        alpha[mask] = 0

        tau_max = np.sum(alpha, axis=0) * dx

        I = I_tau(tau_max, nu, exT)

        F = np.sum(I * dA) / (d**2)
        F = toMilliJansky(F)
        F_Mdot.append(F)

    return massloss, F_Mdot


cell_size = 0.1
high_grid_size = 1000
low_grid_size = 100
high_cells = high_grid_size / cell_size
low_cells = low_grid_size / cell_size


# not needed if there was no previous computation
previous_data = pd.read_csv("Radiostars_massloss_3.csv", delimiter=",")
previous_massloss = previous_data["Massloss-rate_constraint"].values.astype(float)

masslossconstraints = []
for i in range(len(ID)):
    if previous_massloss[i] > 1e5 * Mdot_sun:
        # print("high")
        masslossrate, F_Mdot = Fluxcalc(
            Radii[i],
            Masses[i],
            Frequencies[i] * 1e6,
            Distances[i] * pc,
            Rmax=high_grid_size,
            L_grid=int(high_cells),
        )
        limflux = Flux[i]
        limMdot = np.interp(limflux, F_Mdot, masslossrate)
        masslossconstraints.append(np.interp(Flux[i], F_Mdot, masslossrate))

    else:
        # print("low")
        masslossrate, F_Mdot = Fluxcalc(
            Radii[i],
            Masses[i],
            Frequencies[i] * 1e6,
            Distances[i] * pc,
            Rmax=low_grid_size,
            L_grid=int(low_cells),
        )
        limflux = Flux[i]
        limMdot = np.interp(limflux, F_Mdot, masslossrate)
        masslossconstraints.append(np.interp(Flux[i], F_Mdot, masslossrate))


# add massloss constraints and save as new csv file
data["Massloss-rate_constraint"] = masslossconstraints
pd.DataFrame(data).to_csv("Radiostars_massloss_6.csv", index=False)


# below is for plotting proxima centauri of other stars.
exit()

index = -1  # G 131-26
# index = 180  # proxima centauri


masslossrate, F_Mdot = Fluxcalc(
    0.12,
    0.15,
    Frequencies[index] * 1e6,
    pc * 4.24 / 3.26,
)

limflux = Flux[index]
limit1 = np.full((100), limflux)
limMdot = np.interp(limflux, F_Mdot, masslossrate)
limit2 = np.full((100), limMdot)
lum = []
for i, f in enumerate(F_Mdot):
    lum.append(getLum(f, Distances[i]))

plt.figure()
plt.title(ID[index])
# plt.title("Proxima Centauri")  # for proxima centauri specifically
plt.plot(masslossrate / Mdot_sun, F_Mdot)
plt.axhline(limflux, color="k")
# plt.axvline(limMdot / Mdot_sun, color="k")
# plt.plot(masslossrate / Mdot_sun, limit1, "black")
# plt.plot(limit2 / Mdot_sun, F_Mdot, "black")

ax = plt.gca()
ax.set_ylim([0, 11])

plt.xlabel("Mass-loss rate ($\dot{M}_\odot$)")
plt.xscale("log")
plt.ylabel("Radio Luminosity ($W Hz^{-1}$)")
plt.ylabel("Flux density (mJy)")

plt.show()
