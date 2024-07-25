# read the radiostars catalogue csv file
# consider each quantity in this file
# query Simbad database
# create csv file with desired quantities

import pandas as pd
import numpy as np
from astroquery.simbad import Simbad

# radiostars catalogue
data_path = r"C:\Users\Arnau\Documents\Bachelorproject\Python\Radio.dat"


data = pd.read_csv("SRSC_filtered.csv", delimiter=",")
Simbad_name = data["Simbad_ID"].values

frequency = data["Radio_freq_MHz"].values.astype(float)
fluxdensity = data["Radio_I_flux_int"].values.astype(float)

parallax = data["Archival_parallax"].values.astype(float)
d = 1000 / parallax


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


# query Simbad database
Mks_list = []
Mass = []
Simbad.add_votable_fields("flux(K)")
for qrysimbad in Simbad_name:
    if str(qrysimbad) != "nan":
        try:
            Mks = Simbad.query_object(str(qrysimbad))["FLUX_K"].value[0]
            Mks_list.append(Mks)
            Mass.append(getMass(Mks))
        except TypeError:
            Mks_list.append("Error Name")
            Mass.append("Error Name")
    else:
        Mks_list.append("No Name")
        Mass.append("No Name")


def datatocsv(FolderLocation, dataset, filename, columnnames):
    """Make csv file of given data.
    FolderLocation = location of the folder where the csv file is stored
    dataset = numpy 2d array of the data. First index number is the quantity
    filename = csv file name
    columnsnames = list of names of each column for the header, in order!
    """

    # take transpose as first index number selects the column
    df = pd.DataFrame(data=np.transpose(dataset), columns=columnnames)

    fileLocation = FolderLocation + rf"\{filename}.csv"
    df.to_csv(fileLocation)


# add all desired quantities to a single np 2darray
dataset = np.concatenate(
    ([Simbad_name], [frequency], [fluxdensity], [parallax], [d], [Mks_list], [Mass])
)


# label each collumn in the same order as done in np.concatenate
columns = [
    "Simbad ID",
    "Radio_freq_MHz",
    "Radio_I_flux_int",
    "Archival_parallax",
    "distance",
    "Flux K",
    "Mass",
]


# location to store the csv file
folder_path = r"C:\Users\Arnau\Documents\Bachelorproject\Python\QueryResults"

csvfilename = "Radiostars_Complete"

# create csv
datatocsv(folder_path, dataset, filename=csvfilename, columnnames=columns)
