# Read radio stars catalogue csv
# filter the stars including a Simbad or Gaia ID
# query the Simbad database
# if nothing, then query Gaia database
# return csv file of ID, apparent K-magnitude, parallax and whether they are obtained from Simbad or Gaia


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


def getKs(G, G_BPRP):
    Ks = G - (-0.0981 + 2.098 * G_BPRP - 0.1579 * (G_BPRP**2))
    return Ks


def query(simbad_id, gaia_id):

    Ks = []
    plx = []
    source = []
    ID = []
    Simbad.add_votable_fields("plx", "flux(K)")
    for i, id in enumerate(simbad_id):
        if id != "":  # check if there is a simbad id
            if np.str_(id) != "nan":
                if id != simbad_id[i - 1]:  # check if it is not a duplicate
                    obj = Simbad.query_object(id)
                    try:
                        if not obj["FLUX_K"].mask[0]:  # check if K mag can be retrieved
                            if not obj["PLX_VALUE"].mask[
                                0
                            ]:  # check if parallax can be retrieved
                                Ks.append(obj["FLUX_K"].value[0])
                                plx.append(obj["PLX_VALUE"].value[0])
                                source.append("Simbad")
                                ID.append(id)
                        elif gaia_id[i] != "":
                            if np.str_(gaia_id[i]) != "nan":
                                gaia_number = gaia_id[i].replace("Gaia DR3 ", "")

                                qrygaia = f"SELECT phot_g_mean_mag, bp_rp, parallax FROM gaiadr3.gaia_source WHERE source_id = {gaia_number}"
                                job = Gaia.launch_job_async(qrygaia)
                                row = job.get_results()
                                if row["phot_g_mean_mag"][0] < 13:
                                    Ks.append(
                                        getKs(
                                            row["phot_g_mean_mag"][0], row["bp_rp"][0]
                                        )
                                    )
                                    plx.append(row["parallax"][0])
                                    source.append("Gaia")
                                    ID.append(gaia_id[i])
                    except TypeError:
                        if np.str_(gaia_id[i]) != "nan":
                            gaia_number = gaia_id[i].replace("Gaia DR3 ", "")

                            qrygaia = f"SELECT phot_g_mean_mag, bp_rp, parallax FROM gaiadr3.gaia_source WHERE source_id = {gaia_number}"
                            job = Gaia.launch_job_async(qrygaia)
                            row = job.get_results()
                            if row["phot_g_mean_mag"][0] < 13:
                                Ks.append(
                                    getKs(row["phot_g_mean_mag"][0], row["bp_rp"][0])
                                )
                                plx.append(row["parallax"][0])
                                source.append("Gaia")
                                ID.append(gaia_id[i])
    return Ks, plx, source, ID


Ks, plx, source, ID = query(Simbad_ID, Gaia_ID)

pd.DataFrame(
    {"ID": ID, "K_mag": Ks, "Parallax": plx, "source_obtained": source}
).to_csv("mags_plx_combined_id.csv", index=False)


#########################
# old code below
#######################
exit()


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
