# Compute spectrum of different grid sizes and cell size
# saving each plot as a txt with the title as their label to save repeated computations

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


def Fluxcalc(Rstar, Mstar, Mdot, d, Rmax=1500, L_grid=500):
    """Calculate the flux density of a star using the isothermal wind model for the density profile.
    Rstar = radius of the star in solar radii
    Mstar = mass of the star in solar masses
    Mdot = mass-loss rate in g/s
    d = distance in cm
    Rmax = system size in solar radii
    L_grid = number of points on each side

    Returns list of frequencies and corresponding flux densities.
    """

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

    freq = np.logspace(8, 9, 20) * 1.5  # frequency range 144-1655 MHz

    # temperature fixed to value: (example value exT=2e6 K)
    r, u = isothermal_wind(Mstar, Rstar, exT, Rmax, npts=L_grid)

    F_nu = []
    for nu in freq:

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
        F_nu.append(F)

    return freq, F_nu


# add massloss constraints and save as new csv file
# data["Massloss-rate_constraint"] = masslossconstraints
# pd.DataFrame(data).to_csv("Radiostars_massloss_3.csv", index=False)


def plotgridsizes(gridsizes, cellsize, Mdot):
    plt.figure()
    plt.title(
        ID[-1]
        + " | $\dot{M}$: "
        + str(Mdot / Mdot_sun)
        + " & cell size: "
        + str(cellsize)
    )

    gridsizes = sorted(gridsizes, reverse=True)

    for p in gridsizes:
        npts = p / cellsize
        nu, F = Fluxcalc(
            Radii[-1],
            Masses[-1],
            Mdot,
            Distances[-1] * pc,
            L_grid=int(npts),
            Rmax=p,
        )
        if p == max(gridsizes):
            highest_Fres = F

        perc = percdiff(F, highest_Fres)
        pd.DataFrame({"nu": nu, "F": F, "perc": perc}).to_csv(
            rf"textdata\GS{str(p)}CS{str(round(cellsize, 3))}Mdot{str(Mdot/Mdot_sun)}.txt",
            index=False,
            header=False,
        )

        plt.plot(nu, perc, label=f"Grid size: {p}")
        print(f"{p} done")
    # plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Percentage difference")
    # plt.ylabel("Flux density (mJy)")
    plt.legend(loc="lower right")

    # plt.show()


def plotcellsizes(gridsize, cellsizes, Mdot):
    plt.figure()
    plt.title(
        ID[-1]
        + " | $\dot{M}$: "
        + str(Mdot / Mdot_sun)
        + " & $R_{max}$: "
        + str(gridsize)
    )
    for p in cellsizes:
        npts = gridsize / p
        nu, F = Fluxcalc(
            Radii[-1],
            Masses[-1],
            Mdot,
            Distances[-1] * pc,
            L_grid=int(npts),
            Rmax=gridsize,
        )
        if p == min(cellsizes):
            highest_Fres = F

        perc = percdiff(F, highest_Fres)

        pd.DataFrame({"nu": nu, "F": F, "perc": perc}).to_csv(
            rf"textdata\GS{str(gridsize)}CS{str(round(p,3))}Mdot{str(Mdot/Mdot_sun)}.txt",
            index=False,
            header=False,
        )

        plt.plot(nu, perc, label=f"{round(p, 3)} Rstar / cell")
        print(f"{p} done")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Percentage difference")
    plt.legend(loc="lower right")

    # plt.show()


def percdiff(spectrum, highest_res):
    q = []
    for i in range(len(spectrum)):
        q.append(100 * np.abs(spectrum[i] - highest_res[i]) / highest_res[i])
    return q


# example: last star on the list
low_Mdot = 1e2 * Mdot_sun
high_Mdot = 1e6 * Mdot_sun


grid_sizes = np.logspace(1, 3, 3) * 5
cell_sizes = np.linspace(10, 30, 3)

# ideal for high_Mdot
grid_size = 3000  # Rstar
cell_size = 10  # Rstar / npts

# ideal for low_Mdot
# grid_size = 500
# cell_size = 0.5


manual_grid_sizes = [500, 1000, 1500]
# plotgridsizes(manual_grid_sizes, 0.1, 1e6 * Mdot_sun)


manual_cell_sizes = [0.05, 0.1, 1, 5, 10]
masslosses = [1e3, 1e4, 1e5, 1e6]
for mdot in masslosses:
    plotgridsizes(manual_grid_sizes, 0.1, mdot * Mdot_sun)
# plotgridsizes([500, 1000, 1500], 0.1, 1e6 * Mdot_sun)

plt.show()
