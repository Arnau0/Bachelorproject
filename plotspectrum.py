# read txt files from FreqRadiosstars.py and manipulate plot a desired

import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

Mdot_sun = 1.261e12


def getLum(F, d):
    return 4 * np.pi * F * d**2


def percdiff(spectrum, highest_res):
    q = []
    for i in range(len(spectrum)):
        q.append(100 * np.abs(spectrum[i] - highest_res[i]) / highest_res[i])
    return q


# loc = rf"textdata\GS{str(gridsize)}CS{str(round(p,3))}Mdot{str(Mdot/Mdot_sun)}.txt"
# loc = r"textdata\GS1500CS0.1Mdot1000000.0.txt"

masslosses = [1e4, 1e5, 1e6]
manual_cell_sizes = [0.05, 0.1, 1, 5, 10]
manual_grid_sizes = [50, 100, 150]


def varyCS(cellsizes):
    plt.figure()
    for i, mdot in enumerate(masslosses):

        plt.title("$\dot{M}$: " + str(mdot) + " & $R_{max}$: " + str(1000), fontsize=18)

        for cs in cellsizes:

            loc = rf"textdata\GS{str(1000)}CS{str(round(cs,3))}Mdot{str(mdot)}.txt"
            data = pd.read_csv(
                loc,
                delimiter=",",
                header=None,
            )

            nu = data.iloc[:, 0].astype(float)
            F = data.iloc[:, 1].astype(float)
            if cs == min(cellsizes):
                highest_Fres = F

            perc = percdiff(F, highest_Fres)

            plt.plot(nu, perc, label=f"Cell size: {round(cs, 3)}")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("Frequency (Hz)", fontsize=18)
        # plt.ylabel("Flux density (mJy)", fontsize=18)
        plt.ylabel("Percentage difference", fontsize=18)
        plt.legend(loc="lower right", fontsize=16)
        plt.show()


def varyGS(gridsizes, masslosses):
    plt.figure()
    for i, mdot in enumerate(masslosses):

        plt.title("$\dot{M}$: " + str(mdot) + " & Cell size: " + str(0.1), fontsize=18)

        gridsizes = sorted(gridsizes, reverse=True)
        for gs in gridsizes:

            loc = rf"textdata\GS{str(gs)}CS{str(0.1)}Mdot{str(mdot)}.txt"
            data = pd.read_csv(
                loc,
                delimiter=",",
                header=None,
            )

            nu = data.iloc[:, 0].astype(float)
            F = data.iloc[:, 1].astype(float)
            if gs == max(gridsizes):
                highest_Fres = F

            perc = percdiff(F, highest_Fres)
            plt.plot(nu, perc, label=f"Grid size: {gs}")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("Frequency (Hz)", fontsize=18)
        # plt.ylabel("Flux density (mJy)", fontsize=18)
        plt.ylabel("Percentage difference", fontsize=18)
        plt.legend(loc="lower left", fontsize=16)
        plt.show()


masslosses = [1e3]

# varyCS(manual_cell_sizes)
varyGS(manual_grid_sizes, masslosses)
