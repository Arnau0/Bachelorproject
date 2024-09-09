# Read csv of ID, apparent K-magnitude, parallax and whether they are obtained from Simbad or Gaia
# convert apparent magnitudes to absolute magnitudes
# calculate mass, radius and distance
# return given csv file including mass, radius, distance and absolute K-magnitude
# included are histograms

import pandas as pd
import numpy as np
import math as math
import matplotlib.pyplot as plt


def getMass(M_ks):
    """Calculate the mass of a star using the semi-empirical mass-K magnitude relation from Mann et al. (2015).
    Eq (10)
    M_ks = K magnitude observed from the star
    returns mass in stellar masses
    """
    a = 0.5858
    bX = 0.3872 * M_ks
    cX2 = -0.1217 * M_ks**2
    dX3 = 0.0106 * M_ks**3
    eX4 = -2.7262e-4 * M_ks**4
    return a + bX + cX2 + dX3 + eX4


def getRadius(M_ks):
    """Calculate radius of a star using semi-empirical radius-K magnitude relation from Mann et al.
    Eq (4)
    M_ks = K magnitude observed from the star
    returns mass in stellar radii
    """
    a = 1.9515
    bX = -0.3520 * M_ks
    cX2 = 0.01680 * M_ks**2
    return a + bX + cX2


def getAbsMag(m, d):
    return m - 5 * (np.log10(np.abs(d)) - 1)


data = pd.read_csv(
    r"C:\Users\Arnau\Documents\Bachelorproject\Python\mags_plx_combined_id.csv",
    delimiter=",",
)
data = data[data["Parallax"] != "--"]
data = data[data["K_mag"] != "--"]


Ks = data["K_mag"].values.astype(float)
plx = data["Parallax"].values.astype(float)


MKs_1 = getAbsMag(Ks, 1000 / plx)


plt.figure()
plt.hist(MKs_1, bins=30)
plt.title("Gaia/Simbad query")
plt.xlabel("Absolute K magnitude")
plt.ylabel("Number of stars")
plt.text(-10, 75, f"230 stars within domain", c="black")
plt.axvline(4.6, color="k", linestyle="dashed")
plt.axvline(9.8, color="k", linestyle="dashed")

plt.show()


data = data[(getAbsMag(Ks, 1000 / plx) >= 4.6) * (getAbsMag(Ks, 1000 / plx) < 9.8)]


Ks = data["K_mag"].values.astype(float)
plx = data["Parallax"].values.astype(float)
source = data["source_obtained"].values.astype(str)
ID = data["ID"].values.astype(str)


d = 1000 / plx  # [1000 / x for x in plx]
MKs = getAbsMag(Ks, d)
Mstar = getMass(MKs)
Rstar = getRadius(MKs)


# find F_I and frequency for each ID in radiostars catalogue
data = pd.read_csv(
    r"C:\Users\Arnau\Documents\Bachelorproject\Python\SRSC_filtered.csv",
    delimiter=",",
)

Fstar = []
nustar = []
for i, id in enumerate(ID):
    if source[i] == "Gaia":
        data_source = data[data["GaiaDR3_ID"] == id]
        F_I = min(data_source["Radio_I_flux_int"])
        nu = data_source[data_source["Radio_I_flux_int"] == F_I][
            "Radio_freq_MHz"
        ].values.astype(float)[0]
        Fstar.append(F_I)
        nustar.append(nu)
    else:
        try:
            data_source = data[data["Simbad_ID"] == id]
            F_I = min(data_source["Radio_I_flux_int"])

            nu = data_source[data_source["Radio_I_flux_int"] == F_I][
                "Radio_freq_MHz"
            ].values.astype(float)[0]
            Fstar.append(F_I)
            nustar.append(nu)
        except IndexError:
            Fstar.append("NoMeasurement")
            nustar.append("NoMeasurement")


# create csv if desired
data = pd.DataFrame(
    {
        "ID": ID,
        "Mass": Mstar,
        "Radius": Rstar,
        "Distance": d,
        "Absolute_K_mag": MKs,
        "Source_obtained": source,
        "Flux_I": Fstar,
        "Frequency": nustar,
    }
)
data = data[data["Flux_I"] != "NoMeasurement"]

# get csv file
# data.to_csv("Radiostars_for_flux.csv", index=False)


exit()

# Histograms
plt.figure()
plt.hist(MKs, bins=30)
plt.title("Gaia/Simbad query")
plt.xlabel("Absolute K magnitude")
plt.ylabel("Number of stars")
plt.text(7.7, 14, f"{len(MKs)} stars within domain", c="black")


fig = plt.figure()

fig.add_subplot(1, 2, 1)
plt.hist(Mstar, bins=20)
plt.title("Gaia/Simbad query: Mass")
plt.xlabel("Mass ($M_\odot$)")


fig.add_subplot(1, 2, 2)
plt.hist(Rstar, bins=20)
plt.title("Gaia/Simbad query: Radius")
plt.xlabel("Radius ($R_\odot$)")


plt.figure()
plt.plot(Mstar, Rstar, "bo")
plt.ylabel("Radius ($R_\odot$)")
plt.xlabel("Mass ($M_\odot$)")
plt.show()
