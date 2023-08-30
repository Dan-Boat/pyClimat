# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 13:33:18 2023

@author: dboateng
"""

import pandas as pd 
import numpy as np
import os 
import pygmt


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

#region = [-7, 26, 39, 55]
region = [5, 12, 44.5, 48]
fig = pygmt.Figure()
topo_data = '@earth_relief_30s'

# pygmt.makecpt(
#     cmap='topo',
#     series='-8000/8000/1000',
#     continuous=True
# )

fig.grdimage(grid=topo_data, region=region, projection="M15c", 
             frame=True,)

fig.coast(shorelines=True, borders=["1/1p,white"])       

font = "15p,Helvetica-Bold"
#West/South/East/North coordinates

fig.plot(y=47.012, x=7.955, style="a0.5c", pen="1p,black", fill="darkorange", label="Fontannen")



fig.plot(y=47.285, x=8.978, style="c0.5c", pen="1p,black", fill="darkred", label="Jona")
fig.plot(y=47.233, x=8.627, style="d0.5c", pen="1p,black", fill="seagreen", label="Aabach")
fig.plot(y=46.166, x=8.1466, style="n0.5c", pen="1p,black", fill="lightseagreen", label="SFZ")

fig.plot(x=7.45, y= 46.95, style="c0.2c", pen="1p,black", fill="black")
fig.text(x=7.45, y= 46.95+ 0.15, text="Bern", font=font)

fig.plot(x=8.54, y= 47.38, style="c0.2c", pen="1p,black", fill="black")
fig.text(x=8.54, y= 47.38+ 0.15, text="Zürich", font=font)

# fig.plot(data=np.array([[1,44,8,46.5]]), style='r+s', pen="3p,red")
# fig.plot(data=np.array([[5,46.5,16,50.5]]), style='r+s', pen="3p,black")
# fig.plot(data=np.array([[8,43,16,46.5]]), style='r+s', pen="3p,blue")

#fig.plot(x=[-5, 25], y=[46.25, 46.25], pen="3p,black,-")
#fig.plot(x=[10, 10], y=[40, 54], pen="3p,red,-")




# fig.plot(data=np.array([[-5,45.5,25,46]]), style='r+s', pen="3p,black,-")
# fig.plot(data=np.array([[10,40,10.5,54]]), style='r+s', pen="3p,red,-")


fig.colorbar(frame='+l"Topography"')
fig.legend()
fig.savefig(fname= os.path.join(path_to_plots, "topo_with_proxy_location.pdf"), crop=True, dpi=600)

fig.show()