# plots of the initial results

import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

Mdot_sun = 1.261e12


def getLum(F, d):
    return 4 * np.pi * F * d**2


data = pd.read_csv(
    r"C:\Users\Arnau\Documents\Bachelorproject\Python\Radiostars_massloss_6.csv",
    delimiter=",",
)

ID = data["ID"].values.astype(str)
Masses = data["Mass"].values.astype(float)  # in solar masses
masslossrates = data["Massloss-rate_constraint"].values.astype(float)
distances = data["Distance"].values.astype(float)
flux = data["Flux_I"].values.astype(float)
radii = data["Radius"].values.astype(float)

lum = getLum(flux / int(1e26), distances)


plt.figure()
plt.scatter(Masses, masslossrates / (Mdot_sun), c=flux, norm=LogNorm())
plt.yscale("log")
cbar = plt.colorbar()
# plt.xscale("log")
# plt.ylabel("Mass ($M_\odot$)")
plt.ylabel("Massloss rate ($\dot{M}_\odot$)", fontsize=16)
plt.xlabel("Mass ($M_\odot$)", fontsize=16)
cbar.set_label(label="Flux density (mJy)", fontsize=16)
plt.show()


#######################################333
# plot histogram of stars within 25 pc
# including the one from Winters et al. (2019) for which the values is obtained through inspection

plt.figure()
plt.title("Masses within 25 pc")
# counts = np.array([12, 47, 55, 63, 85, 104, 89, 120, 107, 132, 155, 149])
# masses = np.linspace(0.65, 0.1, len(counts))
# plt.step(masses, counts, label="Winters et al. (2019)", where="mid", color="orange")
counts = np.array([12, 47, 55, 63, 85, 104, 89, 120, 107, 132, 155, 149])
masses = np.linspace(0.675, 0.075, len(counts) + 1)

counts = np.array([149, 155, 132, 107, 120, 89, 104, 85, 63, 55, 47, 12])
masses = np.linspace(0.075, 0.675, len(counts) + 1)
binnumber, massedges = np.histogram(
    Masses[distances < 25], range=(0.075, 0.675), bins=12
)

plt.stairs(counts, masses, label="Winters et al. (2019)")
plt.stairs(binnumber, massedges, label="Sample of radio stars")
# plt.hist(
#   Masses[distances < 25],
#   range=(0.075, 0.675),
#   bins=12,
#   label="Sample of radio stars",
#   color="blue",
# )
plt.xlabel("Mass ($M_\odot$)")
plt.ylabel("Number of stars")
plt.yscale("log")

plt.axis([0.675, 0.075, 0, max(counts) + 10])
plt.legend(fontsize=12)


plt.figure()
plt.title("Masses within 25 pc")
massdiff = masses - massedges
ratio = binnumber / counts

avgratio = np.average(ratio)
plt.stairs(ratio, masses)
plt.axhline(avgratio, color="orange")
print(avgratio)

plt.axis([0.69, 0.03, 0, 0.51])
plt.ylabel("Number ratio")
plt.xlabel("Mass ($M_\odot$)")
plt.show()

exit()
#####################################33


plt.figure()
plt.scatter(
    Masses[distances < 25],
    masslossrates[distances < 25] / Mdot_sun,
    c=lum[distances < 25],
    norm=LogNorm(),
)
plt.yscale("log")
plt.show()


plt.figure()
plt.scatter(distances, lum, c=Masses)
# plt.xlabel("Mass ($M_\odot$)")
plt.xlabel("Distance (pc)")
plt.ylabel("Radio Luminosity (W Hz$^{-1}$)")
plt.yscale("log")
plt.xscale("log")
# plt.colorbar(label="Flux density (mJy)")
plt.colorbar(label="Mass ($M_\odot$)")


plt.figure()
# plt.scatter(Masses, masslossrates / (Mdot_sun * radii**2), c=flux, norm=LogNorm())
plt.scatter(distances, Masses, c=flux, norm=LogNorm())
plt.yscale("log")
plt.colorbar(label="Flux density (mJy)")
# plt.xscale("log")
plt.ylabel("Mass ($M_\odot$)")
# plt.ylabel("Massloss rate ($\dot{M}_\odot$)")
plt.xlabel("Distance (pc)")
# plt.show()


plt.figure()
plt.scatter(Masses, masslossrates / (Mdot_sun * radii ** (-2)), c=flux, norm=LogNorm())
plt.yscale("log")
cbar = plt.colorbar()
# plt.xscale("log")
# plt.ylabel("Mass ($M_\odot$)")
plt.ylabel(
    "Massloss rate per unit surface area ($\dot{M}_\odot$/$R_\odot^2$)", fontsize=16
)
plt.xlabel("Mass ($M_\odot$)", fontsize=16)
cbar.set_label(label="Flux density (mJy)", fontsize=16)
# plt.show()


plt.figure()
plt.scatter(
    distances,
    masslossrates / (Mdot_sun),
    c=flux,
    norm=LogNorm(),
)
plt.yscale("log")
cbar = plt.colorbar()
plt.xscale("log")
plt.ylabel("Massloss rate ($\dot{M}_\odot$)", fontsize=16)
plt.xlabel("Distance (pc)", fontsize=16)
cbar.set_label(label="Flux density (mJy)", fontsize=16)
# plt.show()


plt.figure()
plt.scatter(Masses, lum, c=masslossrates / Mdot_sun, norm=LogNorm())
plt.yscale("log")
plt.colorbar(label="Massloss rate ($\dot{M}_\odot$)")
# plt.xscale("log")
plt.ylabel("Radio Luminosity (W Hz$^{-1}$)")
plt.xlabel("Mass ($M_\odot$)")
# plt.show()

plt.figure()
plt.hist(Masses, bins=10)
plt.xlabel("Mass")
# plt.show()

starswithin = sum(distances < 25)


plt.figure()
plt.scatter(Masses, lum, c=masslossrates / Mdot_sun, norm=LogNorm())
plt.yscale("log")
plt.colorbar(label="Massloss rate ($\dot{M}_\odot$)")
# plt.xscale("log")
plt.ylabel("Radio Luminosity (W Hz$^{-1}$)")
plt.xlabel("Mass ($M_\odot$)")
# plt.show()
