# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 13:33:18 2023

@author: dboateng
"""

import pandas as pd 
import numpy as np
import os 
import pygmt


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/plots"

region = [-7, 22, 42, 55]
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

fig.coast(shorelines=True)       


#West/South/East/North coordinates


# fig.plot(data=np.array([[1,44,8,46.5]]), style='r+s', pen="3p,red")
# fig.plot(data=np.array([[5,46.5,16,50.5]]), style='r+s', pen="3p,black")
# fig.plot(data=np.array([[8,43,16,46.5]]), style='r+s', pen="3p,blue")

fig.plot(data=np.array([[-5,45.5,20,46.5]]), style='r+s', pen="3p,black,-")
fig.plot(data=np.array([[10,43,15,54]]), style='r+s', pen="3p,red,-")


fig.colorbar(frame='+l"Topography"')
fig.savefig(fname= os.path.join(path_to_plots, "profile_sections.pdf"), crop=True, dpi=600)

fig.show()