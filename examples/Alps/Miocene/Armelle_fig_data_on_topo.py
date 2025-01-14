# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 18:21:13 2024

@author: dboateng
"""

import pandas as pd 
import numpy as np
import os 
import pygmt

gnip_path = "D:/Datasets/GNIP_data/world/scratch/station_world_overview_5years.csv" 
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

df = pd.read_csv(gnip_path)

# df = df[(df["lat"] >= 35) & (df["lat"] <= 50)]
# df = df[(df["lon"] >= 6) & (df["lon"] <= 10)]

#region = [-7, 26, 39, 55]
region = [5, 20, 40, 50]
fig = pygmt.Figure()
topo_data = '@earth_relief_30s'

# pygmt.makecpt(
#     cmap='topo',
#     series='-8000/8000/1000',
#     continuous=True
# )

fig.grdimage(grid=topo_data, region=region, projection="M15c", 
             frame=True,)
fig.colorbar(frame='+l"Topography"')

fig.coast(shorelines=True, borders=["1/1p,white"])       

font = "15p,Helvetica-Bold"
# Pen
pygmt.makecpt(
    cmap="spectral",
    series=[-15, 0, 0.5],)
    
fig.plot(
    x=df.lon.values,
    y=df.lat.values,
    color=df.d18op.values,
    cmap=True,
    style="c0.5c",
    pen="black",
)

fig.colorbar(
    cmap=True,
    frame="xa1f0.5+lpen color",
    position="JRM",
)


fig.legend()
fig.savefig(fname= os.path.join(path_to_plots, "obs_on_topo.pdf"), crop=True, dpi=600)

fig.show()