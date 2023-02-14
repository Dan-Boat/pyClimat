# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 12:37:09 2023

@author: dboateng
"""

import pandas as pd 
import os 
import pygmt


path_to_data = "D:/Datasets/GNIP_data/station_info_all.csv"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"

df = pd.read_csv(path_to_data)

# region = [df.lon.min() -1,
#           df.lon.max() + 1,
#           df.lat.min() -1, 
#           df.lat.max() + 1]


#region = [-25, 60, 30, 70]
region = [-10, 20, 35, 55]
fig = pygmt.Figure()

fig.basemap(region=region, projection="M15c", frame=True)
fig.coast(land="black", water="skyblue")
fig.plot(x=df.lon, y=df.lat, style="c0.3c", color="white", pen="black")
fig.savefig(fname= os.path.join(path_to_plots, "stations_alps.png"))

fig.show()