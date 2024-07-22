# read the radiostars catalogue csv file
# consider each quantity in this file
# query Simbad database
# create csv file with desired quantities

import pandas as pd
import numpy as np
import csv
from astroquery.simbad import Simbad

# radiostars catalogue
data_path = r"C:\Users\Arnau\Documents\Bachelorproject\Python\Radio.dat"


data = pd.read_csv("SRSC_filtered.csv", delimiter=",")
Simbad_name = data["Simbad_ID"].values

frequency = data["Radio_freq_MHz"].values.astype(float)
fluxdensity = data["Radio_I_flux_int"].values.astype(float)

parallax = data["Archival_parallax"].values.astype(float)
d = 1000 / parallax

# query Simbad database
#
# in progress...
qrysimbad = "Wolf 1225"
result = Simbad.add_votable_fields("flux(K)")
Mks = Simbad.query_object(qrysimbad)["FLUX_K"].value[0]
print(Mks)


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
dataset = np.concatenate(([Simbad_name], [frequency], [fluxdensity]))


# label each collumn in the same order as done in np.concatenate
columns = ["Simbad ID", "Radio_freq_MHz", "Radio_I_flux_int"]


# location to store the csv file
folder_path = r"C:\Users\Arnau\Documents\Bachelorproject\Python\QueryResults"

csvfilename = "test3"

# create csv
# datatocsv(folder_path, dataset, filename=csvfilename, columnnames=columns)
