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
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var

from PIL import Image
import glob


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots_prec = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots/animate/prec"
path_to_plots_temp = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots/animate/temp"
path_to_plots_d18Op = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots/animate/d18Op"


# define paths 
main_path = "D:/Datasets/Model_output_pst/"
lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")


filename_lterm = "1003_1017_1m_mlterm.nc"
# read long-term means

PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)


# extract values and compute long-term means

def extract_vars_and_analysis(data, wiso=None):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    #d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)

    
    #compute climatologies
    temp2_alt = compute_lterm_mean(data=temp2, time="month")
    prec_alt = compute_lterm_mean(data=prec, time="month")
   
    
    
    return_data = {"temperature":temp2_alt, "precipitation":prec_alt,
                   }
    
    return return_data

def extract_vars_and_analysis_anomalies(data, pi_data):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
   
    
    temp2_pi = extract_var(Dataset=pi_data , varname="temp2", units="°C")
    prec_pi = extract_var(Dataset= pi_data , varname="prec", units="mm/month")
   
    
    #compute climatologies difference
    
    temp2_diff = compute_lterm_diff(data_control=temp2_pi, data_main=temp2, time="month")
    prec_diff = compute_lterm_diff(data_control=prec_pi, data_main=prec, time="month")
    
    
    
    return_data = {"temperature":temp2_diff, "precipitation":prec_diff,}
    
    return return_data
    
# read data 
lgm_data = extract_vars_and_analysis(data=LGM_data)
mh_data = extract_vars_and_analysis(data=MH_data)
mplio_data = extract_vars_and_analysis(data=PLIO_data)
pi_data = extract_vars_and_analysis(data=PI_data)

mh_anomaly = extract_vars_and_analysis_anomalies(data=MH_data, pi_data=PI_data)
lgm_anomaly = extract_vars_and_analysis_anomalies(data=LGM_data, pi_data=PI_data)
mplio_anomaly = extract_vars_and_analysis_anomalies(data=PLIO_data, pi_data=PI_data)

import calendar




def animate_prec_africa(data_alt, label):
    

    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    for j,mname in enumerate(mnames):
        ax = plt.gca()
        ax.remove()
        projection = ccrs.PlateCarree()
        
        apply_style(fontsize=28, style=None, linewidth=2.5) 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,10), subplot_kw={"projection":projection})
        
        plot_annual_mean(variable="Precipitation", data_alt=data_alt.get("precipitation").sel(month=j+1), ax=ax,
                          cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                      levels=22, level_ticks=6, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                        left_labels=False, fig=fig, plot_borders=False, 
                        plot_projection=projection, title= label + " (" + mname + ")", orientation="horizontal",
                        cbar_pos= [0.28, 0.03, 0.55, 0.03], center=False)

        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        plt.savefig(os.path.join(path_to_plots_prec,  mname + "africa.png"), format= "png", bbox_inches="tight", dpi=300)

def animate_prec(i, data_alt, label):
    lon = i

    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    for j,mname in enumerate(mnames):
        ax = plt.gca()
        ax.remove()
        projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
        
        apply_style(fontsize=28, style=None, linewidth=2.5) 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})
        
        plot_annual_mean(variable="Precipitation anomalies", data_alt=data_alt.get("precipitation").sel(month=j+1), ax=ax,
                          cmap=YlGnBu, units="mm/month", vmax=450, vmin=0, 
                        levels=22, level_ticks=7, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                        left_labels=False, fig=fig, plot_borders=False, 
                        plot_projection=projection, title= label + " (" + mname + ")", orientation="horizontal",
                        cbar_pos= [0.28, 0.05, 0.45, 0.02], center=False)

        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
        plt.savefig(os.path.join(path_to_plots_prec,  str(i)+ mname + ".png"), format= "png", bbox_inches="tight", dpi=300)
# def animate_d18Op(i):
#     lon = i

#     mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
#     for j,mname in enumerate(mnames):
#         ax = plt.gca()
#         ax.remove()
#         projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
        
#         apply_style(fontsize=28, style=None, linewidth=2.5) 
#         fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})
        
#         plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=Mio278_data.get("d18Op").sel(month=j+1), ax=ax,
#                       cmap=RdYlBu, units="‰", vmax=2, vmin=-30, 
#                     levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
#                     left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
#                     plot_projection=projection, title="MIO 278 ppm (" + mname + ")", orientation="horizontal",
#                     cbar_pos= [0.28, 0.05, 0.45, 0.02], center=False)
    
        
#         fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
#         plt.tight_layout() 
#         plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
#         plt.savefig(os.path.join(path_to_plots_d18Op,  str(i)+ mname + ".png"), format= "png", bbox_inches="tight", dpi=600)
        
        
def animate_temp(i, data_alt, label):
    lon = i

    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    for j,mname in enumerate(mnames):
        ax = plt.gca()
        ax.remove()
        projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
        
        apply_style(fontsize=28, style=None, linewidth=2.5) 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})

        plot_annual_mean(variable="Temperature anomalies", data_alt=data_alt.get("temperature").sel(month=j+1), ax=ax,
                          cmap=Spectral_r, units="°C", vmax=30, vmin=-10, 
                        levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                        left_labels=False, fig=fig, plot_borders=False, 
                        plot_projection=projection, title=label + " (" + mname + ")", orientation="horizontal",
                        cbar_pos= [0.28, 0.05, 0.45, 0.02], center=False)
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
        plt.savefig(os.path.join(path_to_plots_temp,  str(i)+ mname + ".png"), format= "png", bbox_inches="tight", dpi=300)
        
        
def animate_temp_alt(i, data_alt, label):
    lon = i

   
   
    ax = plt.gca()
    ax.remove()
    projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
    
    apply_style(fontsize=28, style=None, linewidth=2.5) 
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})

    plot_annual_mean(variable="Temperature anomalies", data_alt=data_alt.get("temperature").mean(dim="month"), ax=ax,
                      cmap=Spectral_r, units="°C", vmax=30, vmin=-10, 
                    levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=False, 
                    plot_projection=projection, title=label, orientation="horizontal",
                    cbar_pos= [0.28, 0.05, 0.45, 0.02], center=False)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
    plt.savefig(os.path.join(path_to_plots_temp,  str(i) + ".png"), format= "png", bbox_inches="tight", dpi=300)
    
def animate_temp_anomaly_alt(i, data_alt, label):
    lon = i

   
   
    ax = plt.gca()
    ax.remove()
    projection = ccrs.Orthographic(central_latitude=20, central_longitude=lon)
    
    apply_style(fontsize=28, style=None, linewidth=2.5) 
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})

    plot_annual_mean(variable="Temperature anomalies", data_alt=data_alt.get("temperature").mean(dim="month"), ax=ax,
                      cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                    levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=False, 
                    plot_projection=projection, title=label, orientation="horizontal",
                    cbar_pos= [0.28, 0.05, 0.45, 0.02], center=True)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
    plt.savefig(os.path.join(path_to_plots_temp,  str(i) + ".png"), format= "png", bbox_inches="tight", dpi=300)



def create_gif2(path_to_plots, output_name, resize_factor=0.5, duration=0.1):
    # Output GIF file
    output_gif = output_name + ".gif"

    # Get a list of all PNG files in the directory and sort them by creation time
    png_files = sorted(glob.glob(os.path.join(path_to_plots, "*.png")), key=os.path.getctime)

    # Create a list to store the resized image objects
    resized_images = []

    # Calculate the new size
    original_width, original_height = Image.open(png_files[0]).size
    new_width = int(original_width * resize_factor)
    new_height = int(original_height * resize_factor)

    # Loop through the PNG files, resize, and append them to the resized_images list
    for png_file in png_files:
        img = Image.open(png_file)
        resized_img = img.resize((new_width, new_height), Image.ANTIALIAS)
        resized_images.append(resized_img)

    # Save the list of resized images as a GIF
    resized_images[0].save(
        os.path.join(path_to_plots, output_gif),
        save_all=True,
        append_images=resized_images[1:],
        loop=0,
        duration=duration
    )

    # Remove the original PNG files
    for png_file in png_files:
        os.remove(png_file)
        print(f"Deleted: {png_file}")




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



#animate_prec_africa(pi_data, "Pre-Industrial")
        
for i in np.arange(0, 360, 30):
    #animate_temp_alt(i, pi_data, "Pre-Industrial")
    #animate_temp_anomaly_alt(i, mh_anomaly, "MH-PI")
    #animate_temp_anomaly_alt(i, lgm_anomaly, "LGM-PI")
    animate_temp_anomaly_alt(i, mplio_anomaly, "mPlio-PI")
    
#create_gif2(path_to_plots_temp, "temp_annual_mh_anomaly", resize_factor=0.5, duration=1500)
#create_gif2(path_to_plots_temp, "temp_annual_lgm_anomaly", resize_factor=0.5, duration=1500)
create_gif2(path_to_plots_temp, "temp_annual_plio_anomaly", resize_factor=0.5, duration=1500)
#create_gif2(path_to_plots_temp, "temp_pi", resize_factor=0.5, duration=1500)

#create_gif2(path_to_plots_prec, "prec_africa", resize_factor=0.5, duration=1000)
