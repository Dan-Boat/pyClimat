# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 19:07:26 2023

@author: dboateng
"""
import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from pyClimat.plot_utils import *


path_to_results= "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/MH"

df_nao_ea_to_prec = pd.read_csv(os.path.join(path_to_results, "MH_NAO_EA_to_prec.csv"))
df_nao_ea_to_t2m = pd.read_csv(os.path.join(path_to_results, "MH_NAO_EA_to_t2m.csv"))
df_t2m_prec_to_nao = pd.read_csv(os.path.join(path_to_results, "MH_climate_to_NAO.csv"))
df_prec_ea_to_nao = pd.read_csv(os.path.join(path_to_results, "MH_prec_EA_to_NAO.csv"))
df_t2m_ea_to_nao = pd.read_csv(os.path.join(path_to_results, "MH_t2m_EA_to_NAO.csv"))


apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(20, 15), sharex=True)

df_t2m_prec_to_nao.mean().plot(kind="barh", ax=ax1)
df_prec_ea_to_nao.mean().plot(kind="barh", ax=ax2)
df_t2m_ea_to_nao.mean().plot(kind="barh", ax=ax3)

axes = [ax1, ax2, ax3]
titles = ["(a) Y(t2m, Prec) to X(NAO)", "(b) Y(Prec, EA) to X(NAO)",  
          "(c) Y(t2m, EA) to X(NAO)"]

for i, ax in enumerate(axes):
    ax.set_title(titles[i], fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="left")
    ax.axvline(x=0.1, linestyle="-", color="red")
    ax.axvline(x=0.33, linestyle="-", color="blue")
    ax.axvline(x=0.66, linestyle="-", color="magenta")
    
plt.savefig(os.path.join(path_to_results, "causal_regional_means_to_NAO_MH.svg"), format= "svg", 
            bbox_inches="tight", dpi=300)



apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(20, 15), sharex=True)

df_nao_ea_to_prec.mean().plot(kind="barh", ax=ax1)
df_nao_ea_to_t2m.mean().plot(kind="barh", ax=ax2)


axes = [ax1, ax2]
titles = ["(a) Y(NAO, EA) to X(Prec)", "(b) Y(NAO, EA) to X(t2m)"]

for i, ax in enumerate(axes):
    ax.set_title(titles[i], fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="left")
    ax.axvline(x=0.1, linestyle="-", color="red")
    ax.axvline(x=0.33, linestyle="-", color="blue")
    ax.axvline(x=0.66, linestyle="-", color="magenta")
    
plt.savefig(os.path.join(path_to_results, "causal_regional_means_to_Prec_t2m_MH.svg"), format= "svg", 
            bbox_inches="tight", dpi=300)