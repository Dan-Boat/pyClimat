# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 12:37:09 2023

@author: dboateng
"""

import pandas as pd 
import numpy as np
import os 
import pygmt


path_to_data = "D:/Datasets/GNIP_data/processed_all/station_info_all.csv"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"


df = pd.read_csv(path_to_data)

# region = [df.lon.min() -1,
#           df.lon.max() + 1,
#           df.lat.min() -1, 
#           df.lat.max() + 1]


region = [-35, 45, 30, 70]
#region = [-10, 25, 35, 55]
fig = pygmt.Figure()
topo_data = '@earth_relief_30s'

# pygmt.makecpt(
#     cmap='topo',
#     series='-8000/8000/1000',
#     continuous=True
# )

fig.grdimage(grid=topo_data, region=region, projection="M15c", 
             frame=True,)

#fig.basemap(region=region, projection="M15c", frame=True)
fig.coast(shorelines=True)       
fig.plot(x=df.lon, y=df.lat, style="cc", color="white", pen="black",)

#West/South/East/North coordinates


# fig.plot(data=np.array([[5.5,55,35,50]]), style='r+s', pen="3p,blue")
# fig.plot(data=np.array([[5,48,15,54]]), style='r+s', pen="3p,yellow")
# fig.plot(data=np.array([[5.5,43.5,15,48]]), style='r+s', pen="3p,black")


fig.colorbar(frame='+l"Topography"')
fig.savefig(fname= os.path.join(path_to_plots, "stations_alps_to_be_edited.pdf"), crop=True, dpi=720)

fig.show()