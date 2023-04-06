# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 20:11:52 2023

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"
        
        
# read all the required datsets (for winter)
df_era_pcs = pd.read_csv(PCs_ERA_DJF_path , parse_dates=["time"])
df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

# read only the times related to the index sign
same = (df_echam_pcs["1"]*-1 >0)&(df_echam_pcs["2"]> 0) | (df_echam_pcs["1"]*-1 <0)&(df_echam_pcs["2"]< 0)
same_positive = (df_echam_pcs["1"]*-1 >0)&(df_echam_pcs["2"]> 0) 
op = ((df_echam_pcs["1"]*-1 >0)&(df_echam_pcs["2"]< 0) | (df_echam_pcs["1"]*-1 <0)&(df_echam_pcs["2"]> 0))

df_echam_nao_ea_same = df_echam_pcs[same]
df_echam_nao_ea_same_positive = df_echam_pcs[same_positive]
df_echam_nao_ea_op = df_echam_pcs[op]



nao_same_echam = xr.DataArray(df_echam_nao_ea_same[str(1)] *-1, dims="time", coords={"time": df_echam_nao_ea_same["time"]})
nao_same_echam_positive = xr.DataArray(df_echam_nao_ea_same_positive[str(1)] *-1, dims="time", coords={"time": df_echam_nao_ea_same_positive["time"]})

nao_op_echam = xr.DataArray(df_echam_nao_ea_op[str(1)] *-1, dims="time", coords={"time": df_echam_nao_ea_op["time"]})

def perform_correlation(atmos_index, varname, units=None, sig=0.05, lev=None, lev_units=None, save_sval=False,
                        path_to_save=main_path_to_data, filename=None):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data = extract_var(Dataset=echam_data, varname=varname, units=units, Dataset_wiso=echam_wiso,
                       lev=lev, lev_units=lev_units)
    
    data_season = extract_region(data=data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    data_season["time"] = data_season.time.values.astype("datetime64[ns]")
    
    data_index = data_season.sel(time = data_season.time.isin(atmos_index.time.values))
    
    sval, pval, sig = StatCorr(x=data_index, y=atmos_index, dim="time",
                                                 return_sig=True, sig=0.05)
    
    if save_sval:
        results = xr.merge([sval, sig])
        results.to_netcdf(os.path.join(path_to_save, filename))
        
    return sval, pval, sig


d18op_nao_same_sval, d18op_nao_same_pval, d18op_nao_same_sig = perform_correlation(atmos_index=nao_same_echam, 
                                                                    varname="d18op", units="per mil",
                                                                    filename="d18op_nao_same_scc.nc")

d18op_nao_op_sval, d18op_nao_op_pval, d18op_nao_op_sig = perform_correlation(atmos_index=nao_op_echam, 
                                                                    varname="d18op", units="per mil",
                                                                    filename="d18op_nao_op_scc.nc")



t2m_nao_same_sval, t2m_nao_same_pval, t2m_nao_same_sig = perform_correlation(atmos_index=nao_same_echam, 
                                                                    varname="temp2", units="°C",
                                                                    filename="t2m_nao_scc.nc")
t2m_nao_op_sval, t2m_nao_op_pval, t2m_nao_op_sig = perform_correlation(atmos_index=nao_op_echam, 
                                                                    varname="temp2", units="°C",
                                                                    filename="t2m_ea_scc.nc")

prec_nao_same_sval, prec_nao_same_pval, prec_nao_same_sig = perform_correlation(atmos_index=nao_same_echam, 
                                                                    varname="prec", units="mm/month",
                                                                    filename="prec_nao_scc.nc")
prec_nao_op_sval, prec_nao_op_pval, prec_nao_op_sig = perform_correlation(atmos_index=nao_op_echam, 
                                                                    varname="prec", units="mm/month",
                                                                    filename="prec_ea_scc.nc")

def cal_diff_nao_ea_index(varname, units=None, index_same=nao_same_echam_positive, index_op=nao_op_echam,
                          lev=None, lev_units=None):

    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data = extract_var(Dataset=echam_data, varname=varname, units=units, Dataset_wiso=echam_wiso,
                       lev=lev, lev_units=lev_units)
    
    data_season = extract_region(data=data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    data_season["time"] = data_season.time.values.astype("datetime64[ns]")
    
    data_same = data_season.sel(time = data_season.time.isin(index_same.time.values)).mean(dim="time")
    data_op = data_season.sel(time = data_season.time.isin(index_op.time.values)).mean(dim="time")
    
    data_diff = data_op - data_same
    data_diff = data_diff.sortby("lon")
    
    return data_diff
    
   
    
d18op_diff = cal_diff_nao_ea_index(varname ="d18op", units="per mil")   

temp_diff = cal_diff_nao_ea_index(varname ="temp2", units="°C")

prec_diff = cal_diff_nao_ea_index(varname ="prec", units="mm/month")


from pyClimat.plots import plot_correlation, plot_annual_mean

# read the save data

apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
#apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.PlateCarree()
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                 figsize=(28, 25), subplot_kw={"projection": projection})

plot_correlation(variable="Spearman Coefficients", data=d18op_nao_same_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.05, 0.01, 0.30, 0.02],
                 plot_pvalues=True, pvalue_data=d18op_nao_same_sig, bottom_labels=True, title="[A] NAO-$\delta^{18}$Op (EQ)")

plot_correlation(variable="Spearman Coefficients", data=d18op_nao_op_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax2, fig=fig, plot_pvalues=True, pvalue_data=d18op_nao_op_sig, bottom_labels=True,
                 left_labels=True, title="[B] NAO-$\delta^{18}$Op (OP)")

plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_diff, cmap="PiYG", units="‰",
                 ax=ax3, fig=fig, vmax=2, vmin=-2, levels=22, domain="Europe Wide", level_ticks=11, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[C] $\delta^{18}$Op (OP-EQ)",
                 left_labels=True, bottom_labels=True, label_format="%.1f")


plot_correlation(variable="Spearman Coefficients", data=t2m_nao_same_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax4, fig=fig, plot_pvalues=True, pvalue_data=t2m_nao_same_sig, left_labels=True,
                 title="[D] NAO-t2m (EQ)", bottom_labels=True)

plot_correlation(variable="Spearman Coefficients", data=t2m_nao_op_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax5, fig=fig, plot_pvalues=True, pvalue_data=t2m_nao_op_sig, left_labels=True,
                 title="[E] NAO-t2m (OP)", bottom_labels=True)

plot_annual_mean(variable='Temperature', data_alt=temp_diff, cmap=RdBu, units="°C",
                 ax=ax6, fig=fig, vmax=4, vmin=-4, levels=22, domain="Europe Wide", level_ticks=11, 
                 cbar_pos = [0.60, 0.10, 0.30, 0.02], title="[F] t2m (OP-EQ)",
                 left_labels=True, bottom_labels=True, label_format="%.1f")

plot_correlation(variable="Spearman Coefficients", data=prec_nao_same_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax7, fig=fig, plot_pvalues=True, pvalue_data=prec_nao_same_sig, left_labels=True,
                 title="[G] NAO-prec (EQ)")

plot_correlation(variable="Spearman Coefficients", data=prec_nao_op_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax8, fig=fig, plot_pvalues=True, pvalue_data=prec_nao_op_sig, left_labels=True,
                 title="[H] NAO-prec (OP)")

plot_annual_mean(variable='Precipitation', data_alt=prec_diff, cmap=BrBG, units="mm/month",
                 ax=ax9, fig=fig, vmax=50, vmin=-50, levels=22, domain="Europe Wide", level_ticks=9, 
                 cbar_pos = [0.05, 0.10, 0.30, 0.02], title="[I] t2m (OP-EQ)",
                 left_labels=True, bottom_labels=True)



fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.10, hspace=0.001)
plt.savefig(os.path.join(path_to_plots, "composite_nao_ea.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.show()    
