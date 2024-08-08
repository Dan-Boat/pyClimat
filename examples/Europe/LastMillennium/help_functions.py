# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:25:36 2024

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
from matplotlib.legend_handler import HandlerPathCollection
import matplotlib.patches as patches
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates
import cartopy.feature as cfeature

from pyClimat.plot_utils import *

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"

def lipd2df(
    lipd_dirpath,
    pkl_filepath=None,
    col_str=[
        "paleoData_pages2kID",
        "dataSetName",
        "archiveType",
        "geo_meanElev",
        "geo_meanLat",
        "geo_meanLon",
        "year",
        "yearUnits",
        "paleoData_variableName",
        "paleoData_units",
        "paleoData_values",
        "paleoData_proxy",
        "minYear",
        "maxYear",
    ],
    crit = "paleoData_iso2kPrimaryTimeseries" # or "paleoData_useInGlobalTemperatureAnalysis"
):
    """
    Convert a bunch of PAGES2k LiPD files to a `pandas.DataFrame` to boost data loading.

    If `pkl_filepath` isn't `None`, save the DataFrame as a pikle file.

    Parameters:
    ----------
        lipd_dirpath: str
          Path of the PAGES2k LiPD files
        pkl_filepath: str or None
          Path of the converted pickle file. Default: `None`
        col_str: list of str
          Name of the variables to extract from the LiPD files

    Returns:
    -------
        df: `pandas.DataFrame`
          Converted Pandas DataFrame
    """

    # Save the current working directory for later use, as the LiPD utility will change it in the background
    work_dir = os.getcwd()
    # LiPD utility requries the absolute path
    lipd_dirpath = os.path.abspath(lipd_dirpath)
    # Load LiPD files
    lipds = lipd.readLipd(lipd_dirpath)
    # Extract timeseries from the list of LiDP objects
    ts_list = lipd.extractTs(lipds)
    # Recover the working directory
    os.chdir(work_dir)
    # Create an empty pandas.DataFrame with the number of rows to be the number of the timeseries (PAGES2k records),
    # and the columns to be the variables we'd like to extract
    df_tmp = pd.DataFrame(index=range(len(ts_list)), columns=col_str)
    # Loop over the timeseries and pick those for global temperature analysis
    i = 0
    for ts in ts_list:
        if (
            crit in ts.keys()
            and ts[crit] == "TRUE"
        ):
            for name in col_str:
                try:
                    df_tmp.loc[i, name] = ts[name]
                except:
                    df_tmp.loc[i, name] = np.nan
            i += 1
    # Drop the rows with all NaNs (those not for global temperature analysis)
    df = df_tmp.dropna(how="all")
    # Save the dataframe to a pickle file for later use
    if pkl_filepath:
        save_path = os.path.abspath(pkl_filepath)
        print(f"Saving pickle file at: {save_path}")
        df.to_pickle(save_path)
    return df


class SupressOutputs(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


class HandlerGreyPathCollection(HandlerPathCollection):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        grey_handle = super().create_artists(legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans)
        for artist in grey_handle:
            artist.set_facecolor("grey")
        return grey_handle

# plot
def plot_data_available(iso2k, add_variable=False, cmap="RdBu_r", varname=None,
                        legend=False, not_global=True, save_name=None, title=None,
                        add_patches=False):
    # list of markers and colors for the different archive_type
    markers = ["p", "p", "o", "v", "d", "*", "s", "s", "8", "D", "^"]
    
    archive_markers = {
        'Coral': 'o',
        'GlacierIce': 's',
        'LakeSediment': 'D',
        'MarineSediment': '^',
        'MolluskShells': 'v',
        'Speleothem': 'p',
        'Wood': '*'
        }
    
    colors = [
        np.array([1.0, 0.83984375, 0.0]),
        np.array([0.73828125, 0.71484375, 0.41796875]),
        np.array([1.0, 0.546875, 0.0]),
        np.array([0.41015625, 0.41015625, 0.41015625]),
        np.array([0.52734375, 0.8046875, 0.97916667]),
        np.array([0.0, 0.74609375, 1.0]),
        np.array([0.25390625, 0.41015625, 0.87890625]),
        np.array([0.54296875, 0.26953125, 0.07421875]),
        np.array([1, 0, 0]),
        np.array([1.0, 0.078125, 0.57421875]),
        np.array([0.1953125, 0.80078125, 0.1953125]),
    ]
    
    # create the plot
    
    apply_style(fontsize=25, style=None, linewidth=2) 
    
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    #projection = ccrs.PlateCarree()
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 14), subplot_kw={"projection":  
                                                                            projection})
    
    # add plot title
    if title is None:
        ax.set_title(
            f"Iso2k Network (n={len(iso2k)})")
    else:
        ax.set_title(title)
    
    # set the base map
    # ----------------
    ax.set_global()
    
    # add coast lines
    ax.coastlines(resolution = "50m", linewidth=1.5, color="grey")
    
    # add land fratures using gray color
    ax.add_feature(cfeature.LAND, facecolor="gray", alpha=0.3)
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3)
    
    if not_global:
        minLon = -120
        maxLon = 60
        minLat = 25
        maxLat = 85
        ax.set_extent([minLon, maxLon, minLat, maxLat], ccrs.PlateCarree())
    
    # add gridlines for latitude and longitude
    gl=ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1,
                         edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5)
    gl.top_labels = False                  # labesl at top
    gl.right_labels = False
    gl.xformatter = LongitudeFormatter()     # axis formatter
    gl.yformatter = LatitudeFormatter()
    
    
    
    # plot the different archive types
    # -------------------------------
    
    # extract the name of the different archive types
    archive_types = iso2k.archiveType.unique()
    
    
    counts = []
    
    # plot the archive_type using a forloop
    for i, type_i in enumerate(archive_types):
        df = iso2k[iso2k["archiveType"] == type_i]
        # count the number of appearances of the same archive_type
        count = df["archiveType"].count()
        
        counts.append(count)
        
        # generate the plot
        if add_variable:
            scatter =ax.scatter(
                df["geo_meanLon"],
                df["geo_meanLat"],
                marker=archive_markers[type_i],
                c = df[varname].values,
                cmap=cmap,
                edgecolor="k",
                s=160,
                vmax =2, vmin=-2,
                transform=ccrs.Geodetic(),
                label=f"{type_i} (n = {count})",
            )
            if i == 0:
                cbar = fig.colorbar(scatter, ax=ax, orientation='vertical', pad=0.05, shrink=0.25)
                cbar.set_label('Mean (Z-score)')
        else:
            ax.scatter(
                df["geo_meanLon"],
                df["geo_meanLat"],
                marker=archive_markers[type_i],
                color=colors[i],
                edgecolor="k",
                s=150,
                transform=ccrs.Geodetic(),
                label=f"{type_i} (n = {count})",
            )
        
        
    # add legend to the plot
    if legend:
        if add_variable:
            # find a way to make the legend one color
            legend_handles = []
            import matplotlib.lines as mlines
            for i, type_i in enumerate(archive_types):
                legend_handles.append(
                    mlines.Line2D(
                        [],
                        [],
                        color="w",
                        marker=archive_markers[type_i],
                        markerfacecolor="gray",  # Choose a fixed color for the legend marker
                        markeredgecolor="k",
                        markersize=10,
                        linestyle="None",
                        label=f"{type_i} (n = {counts[i]})"
                    )
                )
            ax.legend(
            handles=legend_handles,
            scatterpoints=1,
            bbox_to_anchor=(0, -0.6),
            loc="lower left",
            ncol=2,
            # fontsize=15,
            )
            
        else:
            ax.legend(
                scatterpoints=1,
                bbox_to_anchor=(0, -0.6),
                loc="lower left",
                ncol=2,
                #fontsize=15,
            )
            
    if add_patches:
        ax.add_patch(patches.Rectangle(xy =(-90, 60), width= 75, height=24, ls= "-",
                                       color= "red", transform = ccrs.PlateCarree(), 
                                       fc="None", lw=2.5,))
        
        ax.add_patch(patches.Rectangle(xy =(-15, 48), width= 50, height=36, ls= "-",
                                       color= "blue", transform = ccrs.PlateCarree(), 
                                       fc="None", lw=2.5,))
        
        ax.add_patch(patches.Rectangle(xy =(-15, 30), width= 50, height=18, ls= "-",
                                       color= "magenta", transform = ccrs.PlateCarree(), 
                                       fc="None", lw=2.5,))
        
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)
    if save_name is not None:
        plt.savefig(os.path.join(path_to_plots, save_name), format= "png", bbox_inches="tight")
    

def convert_to_float_array(values):
    return np.array(values, dtype=float)


def mean_zscore_950_1250(years, zscores):
    # Filter indices for the years between 950 and 1250
    indices = (years >= 950) & (years <= 1250)
    if np.any(indices):
        filtered_zscores = zscores[indices]
        return np.mean(filtered_zscores)
    else:
        return np.nan
    
    
def mean_zscore_1650_1850(years, zscores):
    # Filter indices for the years between 950 and 1250
    indices = (years >= 1650) & (years <= 1850)
    if np.any(indices):
        filtered_zscores = zscores[indices]
        return np.mean(filtered_zscores)
    else:
        return np.nan
    
    
def mean_zscore_post_1850(years, zscores):
    # Filter indices for the years between 950 and 1250
    indices = (years >= 1851) 
    if np.any(indices):
        filtered_zscores = zscores[indices]
        return np.mean(filtered_zscores)
    else:
        return np.nan
    
    
# Define a function to compute z-scores
def compute_zscore(values):
    if len(values) > 1:  # To avoid computing z-scores for single value arrays
        return zscore(values)
    else:
        return np.zeros_like(values)  # Return an array of zeros if only one value