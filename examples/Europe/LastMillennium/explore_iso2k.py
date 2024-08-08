# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:52:55 2024

@author: dboateng
"""
import os
from io import StringIO
import sys
import numpy as np
from scipy.stats import zscore
import pandas as pd
import lipd
import pyleoclim as pyleo
import tempfile
import pooch

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates
import cartopy.feature as cfeature
import seaborn as sns

from pyClimat.plot_utils import *

path_to_lipd_iso2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/iso2k1_0_1"
path_to_pickle_iso2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/iso2k1_0_1.pkl"
path_to_lipd_pages2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/pages2k"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"


# with SupressOutputs():
#     iso2k = lipd2df(lipd_dirpath=path_to_lipd_iso2k, pkl_filepath=path_to_pickle_iso2k)

from help_functions import *

df = pd.read_pickle(path_to_pickle_iso2k)

df['paleoData_values'] = df['paleoData_values'].apply(convert_to_float_array)
df['year'] = df['year'].apply(convert_to_float_array)


# Apply the function to compute z-scores
df['paleoData_zscores'] = df['paleoData_values'].apply(compute_zscore)

df['mean_zscore_950_1250'] = df.apply(lambda row: mean_zscore_950_1250(row['year'], row['paleoData_zscores']), axis=1)

df['mean_zscore_1650_1850'] = df.apply(lambda row: mean_zscore_1650_1850(row['year'], row['paleoData_zscores']), axis=1)

df['mean_zscore_post_1850'] = df.apply(lambda row: mean_zscore_post_1850(row['year'], row['paleoData_zscores']), axis=1)



iso2k_NH = df[(df['geo_meanLat'] >= 25) & (df['geo_meanLat'] <= 85) & 
                  (df['geo_meanLon'] >= -120) & (df['geo_meanLon'] <= 50)]

#plot for availability

def plot_zscores():

    
    plot_data_available(iso2k=iso2k_NH, add_variable=False, legend=True, add_patches=True,save_name="data_available.png")
    
    
    iso2k_mca = iso2k_NH.dropna(subset=["mean_zscore_950_1250"])
    plot_data_available(iso2k=iso2k_mca, add_variable=True, varname="mean_zscore_950_1250", 
                        save_name="z_score_mca.png", title="MCA (0950-1250)",legend=True)
    
    
    iso2k_lia = iso2k_NH.dropna(subset=["mean_zscore_1650_1850"])
    plot_data_available(iso2k=iso2k_lia, add_variable=True, varname="mean_zscore_1650_1850", 
                        save_name="z_score_lia.png", title="LIA (1650-1850)",legend=True)
    
    
    iso2k_pi = iso2k_NH.dropna(subset=["mean_zscore_post_1850"])
    plot_data_available(iso2k=iso2k_pi, add_variable=True, varname="mean_zscore_post_1850", 
                        save_name="z_score_pi.png", title="PI (post 1850)", legend=True)
    
    

def extract_region(data, varname="mean_zscore_950_1250"):
    greenland = data[(data['geo_meanLat'] >= 60) & (data['geo_meanLat'] <= 85) & 
                      (data['geo_meanLon'] >= -90) & (data['geo_meanLon'] <= -15)][varname].values
    
    
    south_eu = data[(data['geo_meanLat'] >= 30) & (data['geo_meanLat'] <= 50) & 
                      (data['geo_meanLon'] >= -15) & (data['geo_meanLon'] <= 35)][varname].values
    
    north_eu = data[(data['geo_meanLat'] >= 48) & (data['geo_meanLat'] <= 85) & 
                      (data['geo_meanLon'] >= -15) & (data['geo_meanLon'] <= 35)][varname].values
    
    return [greenland, north_eu, south_eu]


regions = ["Greenland", "North EU", "South EU"]

iso2k_mca = iso2k_NH.dropna(subset=["mean_zscore_950_1250"])

mca = extract_region(data=iso2k_mca, varname="mean_zscore_950_1250")

max_len = max(len(m) for m in mca)

df_mca = pd.DataFrame(index=np.arange(max_len), columns=regions).assign(Paleoclimate="MCA")

for i,r in enumerate(regions):
    current_len = len(mca[i])
    df_mca.loc[:current_len - 1, r] = mca[i]
    
    
iso2k_lia = iso2k_NH.dropna(subset=["mean_zscore_1650_1850"])

lia = extract_region(data=iso2k_lia, varname="mean_zscore_1650_1850")

max_len = max(len(m) for m in lia)

df_lia = pd.DataFrame(index=np.arange(max_len), columns=regions).assign(Paleoclimate="LIA")

for i,r in enumerate(regions):
    current_len = len(lia[i])
    df_lia.loc[:current_len - 1, r] = lia[i]
    
    
iso2k_pi = iso2k_NH.dropna(subset=["mean_zscore_post_1850"])

pi = extract_region(data=iso2k_pi, varname="mean_zscore_post_1850")

max_len = max(len(m) for m in pi)

df_pi = pd.DataFrame(index=np.arange(max_len), columns=regions).assign(Paleoclimate="PI")

for i,r in enumerate(regions):
    current_len = len(pi[i])
    df_pi.loc[:current_len - 1, r] = pi[i]
    

cdf = pd.concat([df_mca, df_lia, df_pi])
mdf = pd.melt(cdf, id_vars=["Paleoclimate"],var_name="zscore")

mdf_mean = mdf.groupby(["Paleoclimate", "zscore"]).mean()


def plot_distribution(varname, units, data, hue, ax=None, path_to_plots=None, filename=None,
                       colors=None, xlabel=True, ylabel=True, title=None, ax_legend=True,
                       ymin=None, ymax=None, points_data=None,):
    
    regions = ["Greenland", "North EU", "South EU"]
    
    apply_style(fontsize=25, style=None, linewidth=3)
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
        
    boxplot = sns.violinplot(data=data, x="Paleoclimate", y="value", saturation=0.8, ax=ax,
                          hue=hue, density_norm="width", fill=False, linewidth=3) # count, width, area
    
    if points_data is not None:
        
        pointplot = sns.pointplot(data=points_data, x="Paleoclimate", y="value", ax=ax,
                              hue=hue, markers="x", linestyles="", scale=2,
                              dodge=True, legend=False,)
    
    if colors is not None:
        
        for patch, color in zip(boxplot["boxes"], colors):
            patch.set_facecolor(color)
            
            
    if ylabel:
        ax.set_ylabel(varname + " [" + units + "]", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel is not None:
        ax.set_xlabel("Paleoclimate", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=24,
                  bbox_to_anchor=(0.01, 1.05, 1, 0.102,), loc=3, borderaxespad=0,
                  ncol=4)
    else:
        ax.legend([],[], frameon=False)
        
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 28, "fontweight":"bold"}, loc="center")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
    if path_to_plots is not None:
        plt.savefig(os.path.join(path_to_plots, filename), bbox_inches="tight", format= "png")

apply_style(fontsize=25, style=None, linewidth=3)

plot_distribution(varname="Z-score", units=".", data=mdf, hue="zscore", ymax=3, ymin=-3,
                  ax_legend=True,xlabel=True, title="Mean Z-score (Ref to 095 to 2000)",
                  path_to_plots=path_to_plots, filename="regions.png", points_data=mdf_mean,
                  )
# one_data = iso2k_NH.iloc[66]

# test = pyleo.Series(time=one_data["year"], value=one_data["paleoData_zscores"])
# test.plot()

#plot_data_available(iso2k=df, add_variable=True, varname="mean_zscore_950_1250")
    
    
