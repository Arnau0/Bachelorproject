import pandas as pd
import numpy as np
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad

data = pd.read_csv("SRSC_filtered.csv", delimiter=",")
Gaia_ID = data["GaiaDR3_ID"].values.astype(str)
# Gaia_ID = [x.replace("Gaia DR3 ", "") for x in Gaia_ID]
Simbad_ID = data["Simbad_ID"].values.astype(str)
plx_full = data["Archival_parallax"].values

ID = []
plx = []

# making the list unique for simbad & gaia id's combined
# knowing that the last id in the csv is a duplicate of the second to last id!
for i in range(len(Simbad_ID) - 1):
    if str(Simbad_ID[i]) == "":
        if str(Gaia_ID[i]) != "":
            if Gaia_ID[i] != Gaia_ID[i + 1]:
                ID.append(Gaia_ID[i])
                plx.append(plx_full[i])

    elif str(Simbad_ID[i]) == "nan":
        if str(Gaia_ID[i]) != "":
            if Gaia_ID[i] != Gaia_ID[i + 1]:
                ID.append(Gaia_ID[i])
                plx.append(plx_full[i])

    elif str(Simbad_ID[i]) == "V* BQ Hyi A":
        # specific star in a binary (or more) that isn't found in the simbad database
        # save the Gaia DR3 id instead
        ID.append(Gaia_ID[i])
        plx.append(plx_full[i])

    else:
        if Simbad_ID[i] != Simbad_ID[i + 1]:
            ID.append(Simbad_ID[i])
            plx.append(plx_full[i])

print(np.shape(ID))


def getKs(G, G_BPRP):
    Ks = G - (-0.0981 + 2.098 * G_BPRP - 0.1579 * (G_BPRP**2))
    return Ks


def query(source_ids):

    Ks = []

    Simbad.add_votable_fields("flux(K)")
    for id in source_ids:
        if "Gaia DR3 " in str(id):
            number_id = id.replace("Gaia DR3 ", "")
            qrygaia = f"SELECT phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source WHERE source_id = {number_id}"
            job = Gaia.launch_job_async(qrygaia)
            row = job.get_results()
            Ks_value = getKs(row["phot_g_mean_mag"][0], row["bp_rp"][0])
            Ks.append(Ks_value)
        else:
            Ks_value = Simbad.query_object(str(id))["FLUX_K"].value[0]
            Ks.append(Ks_value)
    return Ks


Ks = query(ID)
print(Ks)
pd.DataFrame({"ID": ID, "Ks": Ks, "Parallax": plx}).to_csv(
    "mags_combined_id.csv", index=False
)


def queryGaia(source_ids):
    """Query the given source ID's and return a pd dataframe of corresponding G and G_bp-Grp magnitudes.
    source_id = list of Gaia DR3 ID's
    """

    ID = []
    G = []
    bprp = []
    for source_ID in source_ids:
        qrygaia = f"SELECT source_id, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source WHERE source_id = {source_ID}"
        job = Gaia.launch_job_async(qrygaia)
        row = job.get_results()
        ID.append(str(row["source_id"][0]))
        G.append(float(row["phot_g_mean_mag"][0]))
        bprp.append(float(row["bp_rp"]))

    dataset = np.concatenate(([ID], [G], [bprp]))
    columnnames = ["source_id", "phot_g_mean_mag", "bp_rp"]
    return dataset, columnnames


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


# location to store the csv file
folder_path = r"C:\Users\Arnau\Documents\Bachelorproject\Python\QueryResults"

csvfilename = "Radiostars_GaiaSimbad"  # without ".csv"

# dataset, columns = queryGaia(GaiaDR3_IDs)

# create csv
# datatocsv(folder_path, dataset, filename=csvfilename, columnnames=columns)
