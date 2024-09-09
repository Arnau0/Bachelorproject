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


d = [1000 / x for x in plx]
MKs = getAbsMag(Ks, d)
MKs = MKs[(MKs >= 4.5) * (MKs < 9.5)]
Mstar = getMass(MKs)
Rstar = getRadius(MKs)


plt.figure()
plt.hist(MKs, bins=30)
plt.title("Gaia/Simbad query")
plt.xlabel("Absolute K magnitude")
plt.text(7.7, 14, f"{len(MKs)} stars within domain", c="black")


fig = plt.figure()
# fig.title("Gaia/Simbad query")

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
