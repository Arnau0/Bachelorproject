import pandas as pd
import numpy as np
import math as math
import matplotlib.pyplot as plt


def getMass(M_ks):
    """Calculate the mass of a star using the semi-empirical mass-K magnitude relation from Mann et al. (2015).
    M_ks = K magnitude observed from the star
    returns mass in stellar radii
    """
    a = 0.5858
    bX = 0.3872 * M_ks
    cX2 = -0.1217 * M_ks**2
    dX3 = 0.0106 * M_ks**3
    eX4 = -2.7262e-4 * M_ks**4
    return a + bX + cX2 + dX3 + eX4


def getAbsMag(m, d):
    return m - 5 * (np.log10(np.abs(d)) - 1)


data = pd.read_csv(
    r"C:\Users\Arnau\Documents\Bachelorproject\Python\mags_combined_id.csv",
    delimiter=",",
)


Ks = []
d = []
for i in range(len(data["ID"])):
    if data["Ks"][i] != "--":
        d_ = 1000 / data["Parallax"][i]
        if str(d_) != "nan":
            Ks.append(float(data["Ks"][i]))
            d.append(1000 / float(data["Parallax"][i]))

print(np.shape(d))
MKs = getAbsMag(Ks, d)
Mstar = getMass(MKs)

plt.figure()
plt.hist(MKs, bins=30)
plt.title("Gaia/Simbad query")
plt.xlabel("Absolute K magnitude")

plt.figure()
plt.hist(Mstar, bins=30)
plt.title("Gaia/Simbad query")
plt.xlabel("Mass ($M_\odot$)")
plt.show()
