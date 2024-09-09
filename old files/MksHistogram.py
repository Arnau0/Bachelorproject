# make histogram of the radiostars with given K magnitude
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv(
    r"C:\Users\Arnau\Documents\Bachelorproject\Python\QueryResults\Radiostars_Complete.csv",
    delimiter=",",
)
Mks_list = data["Flux K"].values

# filter out all strings
Mks_list = Mks_list[Mks_list != "No Name"]
Mks_list = Mks_list[Mks_list != "Error Name"]
Mks_list = [float(x) for x in Mks_list]


# plot histogram
plt.hist(Mks_list, bins=30)
plt.xlabel("K magnitude")
plt.show()
