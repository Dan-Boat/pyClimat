# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 17:16:08 2023

@author: dboateng
"""
import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
import matplotlib.animation as animation
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var

from PIL import Image
import glob


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots_prec = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/animate/prec"
path_to_plots_temp = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/animate/temp"
path_to_plots_d18Op = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/animate/d18Op"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)


# extract values and compute long-term means

def extract_vars_and_analysis(data, wiso, pi_data, pi_wiso):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    
    temp2_pi = extract_var(Dataset=pi_data , varname="temp2", units="°C")
    prec_pi = extract_var(Dataset= pi_data , varname="prec", units="mm/month")
    d18op_pi = extract_var(Dataset=pi_data , varname="d18op", units="per mil", Dataset_wiso= pi_wiso)
    
    #compute climatologies difference
    
    temp2_diff = compute_lterm_diff(data_control=temp2_pi, data_main=temp2, time="month")
    prec_diff = compute_lterm_diff(data_control=prec_pi, data_main=prec, time="month")
    d18op_diff = compute_lterm_diff(data_control=d18op_pi, data_main=d18op, time="month")
    
    
    return_data = {"temperature":temp2_diff, "precipitation":prec_diff, "d18Op":d18op_diff,}
    
    return return_data
    
# read data 
Mio278_data = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso, pi_data=W1E1_PI_data, 
                                        pi_wiso=W1E1_PI_wiso)
Mio450_data = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso, pi_data=W1E1_PI_data, 
                                        pi_wiso=W1E1_PI_wiso)

import calendar




def animate_prec(i):
    lon = i

    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    for j,mname in enumerate(mnames):
        ax = plt.gca()
        ax.remove()
        projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
        
        apply_style(fontsize=28, style=None, linewidth=2.5) 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})
        
        # plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=Mio278_data.get("d18Op").sel(month=j+1), ax=ax,
        #              cmap=GryBr_r, units="‰", vmax=10, vmin=-10, 
        #             levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
        #             left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
        #             plot_projection=projection, title="MIO 278ppm - PI (" + mname + ")", orientation="horizontal",
        #             cbar_pos= [0.28, 0.05, 0.45, 0.02])
        
        plot_annual_mean(variable="Precipitation anomalies", data_alt=Mio450_data.get("precipitation").sel(month=j+1), ax=ax,
                          cmap=PrecAno, units="mm/month", vmax=300, vmin=-300, 
                        levels=22, level_ticks=7, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                        left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                        plot_projection=projection, title="MIO 450ppm - PI (" + mname + ")", orientation="horizontal",
                        cbar_pos= [0.28, 0.05, 0.45, 0.02])

        # plot_annual_mean(variable="Temperature anomalies", data_alt=Mio278_data.get("temperature").sel(month=j+1), ax=ax,
        #                  cmap="RdBu_r", units="°C", vmax=15, vmin=-15, 
        #                 levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
        #                 left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
        #                 plot_projection=projection, title="MIO 278ppm - PI (" + mname + ")", orientation="horizontal",
        #                 cbar_pos= [0.28, 0.05, 0.45, 0.02])
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
        plt.savefig(os.path.join(path_to_plots_prec,  str(i)+ mname + ".png"), format= "png", bbox_inches="tight", dpi=600)


def animate_d18Op(i):
    lon = i

    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    for j,mname in enumerate(mnames):
        ax = plt.gca()
        ax.remove()
        projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
        
        apply_style(fontsize=28, style=None, linewidth=2.5) 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})
        
        plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=Mio450_data.get("d18Op").sel(month=j+1), ax=ax,
                      cmap=GryBr_r, units="‰", vmax=10, vmin=-10, 
                    levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                    plot_projection=projection, title="MIO 450ppm - PI (" + mname + ")", orientation="horizontal",
                    cbar_pos= [0.28, 0.05, 0.45, 0.02])
    
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
        plt.savefig(os.path.join(path_to_plots_d18Op,  str(i)+ mname + ".png"), format= "png", bbox_inches="tight", dpi=600)
        
        
def animate_temp(i):
    lon = i

    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    for j,mname in enumerate(mnames):
        ax = plt.gca()
        ax.remove()
        projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
        
        apply_style(fontsize=28, style=None, linewidth=2.5) 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})

        plot_annual_mean(variable="Temperature anomalies", data_alt=Mio450_data.get("temperature").sel(month=j+1), ax=ax,
                          cmap="RdBu_r", units="°C", vmax=30, vmin=-30, 
                        levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                        left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                        plot_projection=projection, title="MIO 450ppm - PI (" + mname + ")", orientation="horizontal",
                        cbar_pos= [0.28, 0.05, 0.45, 0.02])
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
        plt.savefig(os.path.join(path_to_plots_temp,  str(i)+ mname + ".png"), format= "png", bbox_inches="tight", dpi=600)



def creat_git(path_to_plots, output_name):
    
    # Directory containing PNG files
    output_gif = output_name + ".gif"
    
    # Get a list of all PNG files in the directory
    #png_files = glob.glob(os.path.join(path_to_plots, "*.png"))
    
    # Get a list of all PNG files in the directory and sort them by creation time
    png_files = sorted(glob.glob(os.path.join(path_to_plots, "*.png")), key=os.path.getctime)
    
    # Sort the list of PNG files by filename, if necessary
    #png_files.sort()
    
    # Create a list to store the image objects
    images = []
    
    # Loop through the PNG files and append them to the images list
    for png_file in png_files:
        img = Image.open(png_file)
        images.append(img)
    
    # Save the list of images as a GIF
    images[0].save(os.path.join(path_to_plots,output_gif), save_all=True, 
                   append_images=images[1:], loop=0, duration=0.1)
    
    for png_file in png_files:
        os.remove(png_file)
        print(f"Deleted: {png_file}")
        
for i in np.arange(0, 360, 60):
    animate_prec(i)
    animate_d18Op(i)
    animate_temp(i)
    
    
creat_git(path_to_plots_temp, "MIO_450_temp")
creat_git(path_to_plots_prec, "MIO_450_prec")
creat_git(path_to_plots_d18Op, "MIO_450_d18Op")