# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 22:27:38 2023

@author: dboateng
"""

# importing packages 
import xarray as xr
import os
import pandas as pd
import numpy as np
from scipy import stats
from eofs.xarray import Eof

#importing other routines

try:
    from .data import *
    from .utils import vert_coord_convertion
except:
    from data import *
    from utils import vert_coord_convertion
    
    
    

def extract_var(Dataset, varname, units=None, Dataset_wiso=None, other_data=None, lev_units=None, lev=None,
                use_PDB=False):
    """
    This function extracts some defined variables from a netCDF file. Moreover, if the variable require calculation 
    or unit conversion, the user speficication can be pass to such task. For example, echam out put the differrent 
    component of precipitation (convective, large scale, etc). Therefore extraction of total precipitation would require
    the calculation of the total precipitation. Example of the defined variables: 

    * "temp2" --  near surface temperature
    * "prec" -- total precipitation 
    * "d18op" -- O^18 isotopic composition in precipitation 
    * "d18ov" -- O^18 isotopic composition in vapour
    * "relhum" -- relative humidity
    * "elev"  -- topography or elevation
    "slm" -- mean sea level pressure
    * "evap" -- evaporation
    * "u, v, omega" -- zonal, meridoinal, and vertical velocity
    
    and others 

    example
    
    data = xr.open_dataset(path_to_data)
    temp = extract_var(data, varname= "temp2", unit= "°C",) ---> extract the t2m variable and convert it from K to °C


    Parameters
    ----------
    Dataset : TYPE: Dataset
        DESCRIPTION. Processed output containing all variables or certain ones
    varname : TYPE: STR
        DESCRIPTION: Name of variable (eg. temp2:Tempeature, prec:Total precipitation, d18op: d18o in precipitation)
    units : TYPE, optional 
        DESCRIPTION. The default is None. Units of the variable for eg. °C for temperature or mm/month for precipitation
    Dataset_wiso : TYPE, optional
        DESCRIPTION. The default is None. Wiso Dataset is required for calcullating d18op and d18ov   
    other_data: TYPE: datarray, optional 
        DESCRIPTION: Additional data required for calculating units eg. Pa/s--> m/s require temperature data on pressure levels (in K)
    lev_units: TYPE: str, optional 
        DESCRIPTION: Units of vertical levels eg. hPa
    lev: TYPE: int, optional 
        DESCRIPTION: The vertical height required for the variable, eg. 500 hPa for winds

    Returns
    -------
    var_data : TYPE: datarray
        DESCRIPTION: returns the variable for etracted or calculated from the Dataset (notmally monthly long-term climatologies)

    """
    
    #Temperature
    if hasattr(Dataset, "mlev"):
        Dataset = Dataset.rename({"mlev":"lev"})
        
    if Dataset_wiso is not None:
        if hasattr(Dataset, "mlev"):
            
            Dataset_wiso = Dataset_wiso.rename({"mlev":"lev"})
        
    if varname == "temp2" or  varname == "sst":
        var_data = Dataset[varname]
        if units is not None:
            if units == "°C":
                var_data = var_data - 273.15 # convert temperature to dec C
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kelvin")
                
    elif varname == "tsurf":
        var_data = Dataset[varname]
        if units is not None:
            if units == "°C":
                var_data = var_data - 273.15 # convert temperature to dec C
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kelvin")
                
    elif varname == "tsoil":
        var_data = Dataset[varname]
        var_data = var_data.sel(belowsurface=698)
        
        if units is not None:
            if units == "°C":
                var_data = var_data - 273.15 # convert temperature to dec C
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kelvin")
        
    elif varname == "ws":
        print(".....Extracting soil wetness......")
        var_data = Dataset[varname]
       
        
    #Total precipitation amount
    elif varname == "prec":
        var_data = Dataset["aprl"] + Dataset["aprc"]
        if units is not None:
            if units == "mm/month":
                var_data = var_data *60*60*24*30 #mm/month
                var_data.attrs["units"] = units
            elif units == "mm/a":
                var_data = var_data *60*60*24*365 #mm/a
                var_data.attrs["units"] = units
                
            elif units == "mm/day":
                var_data = var_data *60*60*24 #mm/day
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kg/m²s")
                
                
    #Convective precipitation amount
    elif varname == "aprc":
        var_data = Dataset["aprc"]
        if units is not None:
            if units == "mm/month":
                var_data = var_data *60*60*24*30 #mm/month
                var_data.attrs["units"] = units
            elif units == "mm/a":
                var_data = var_data *60*60*24*365 #mm/a
                var_data.attrs["units"] = units
                
            elif units == "mm/day":
                var_data = var_data *60*60*24 #mm/day
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kg/m²s")
                
    #large-scale precipitation amount
    elif varname == "aprl":
        var_data = Dataset["aprl"]
        if units is not None:
            if units == "mm/month":
                var_data = var_data *60*60*24*30 #mm/month
                var_data.attrs["units"] = units
            elif units == "mm/a":
                var_data = var_data *60*60*24*365 #mm/a
                var_data.attrs["units"] = units
                
            elif units == "mm/day":
                var_data = var_data *60*60*24 #mm/day
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kg/m²s")
    
    # d18o in Precipitation
    elif varname == "d18op":
        var_data_echam = Dataset["aprl"] + Dataset["aprc"]
        if Dataset_wiso is not None:
             var_data_wiso = Dataset_wiso["wisoaprl"][:,1,:,:] + Dataset_wiso["wisoaprc"][:,1,:,:]
        SMOW = 0.2228  #20./18*2005.2e-4
        PDB = 0.22969 #
        wiso_wtgs = var_data_wiso / var_data_echam
        var_data = ((wiso_wtgs/SMOW) -1)
        
        if units is not None:
            if units == "per mil":
                var_data = var_data *1000 # convert to permil
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is isotopic ratio")
                
    elif varname == "d18o_carbonate":
        
        # cal d18Ow
        var_data_echam = Dataset["aprl"] + Dataset["aprc"]
        if Dataset_wiso is not None:
             var_data_wiso = Dataset_wiso["wisoaprl"][:,1,:,:] + Dataset_wiso["wisoaprc"][:,1,:,:]
        
        if use_PDB:
            PDB = 0.22969 #
            wiso_wtgs = var_data_wiso / var_data_echam
            var_data = ((wiso_wtgs/PDB) -1)
        else:
            SMOW = 0.2228  #20./18*2005.2e-4
            wiso_wtgs = var_data_wiso / var_data_echam
            var_data = ((wiso_wtgs/SMOW) -1)
        
        # cal ln(alpha) * 1e-3
        soil_temp = Dataset["tsoil"].sel(belowsurface=698)   
        factor = ((18.03*1000)/soil_temp) - 32.23    
        
        
        if units is not None:
            if units == "per mil":
                var_data = var_data *1000 # convert to permil
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is isotopic ratio")
        
        if use_PDB:
            var_data = var_data + factor
            
        else:
            var_data = ((var_data + factor) -30.91)/1.03091
            
                
                
    # dD in Precipitation
    elif varname == "dDp":
        var_data_echam = Dataset["aprl"] + Dataset["aprc"]
        if Dataset_wiso is not None:
             var_data_wiso = Dataset_wiso["wisoaprl"][:,3,:,:] + Dataset_wiso["wisoaprc"][:,3,:,:]
        SMOW = 0.3288  #19./18.*2.*155.76e-3
        wiso_wtgs = var_data_wiso / var_data_echam
        var_data = ((wiso_wtgs/SMOW) -1)
        
        if units is not None:
            if units == "per mil":
                var_data = var_data *1000 # convert to permil
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is isotopic ratio")
                
    # d18o in water vapour
    elif varname == "d18ov":
        if lev is not None:
            var_data_echam = Dataset["q"].sel(lev=lev)
            if Dataset_wiso is not None:
                var_data_wiso = Dataset_wiso["q18o"].sel(lev=lev)
        else:
            var_data_echam = Dataset["q"][:,30,:,:]
            if Dataset_wiso is not None:
                 var_data_wiso = Dataset_wiso["q18o"][:,30,:,:]
                 
        SMOW = 0.2228
        wiso_wtgs = var_data_wiso / var_data_echam
        var_data = ((wiso_wtgs/SMOW) -1)
        
        if units is not None:
            if units == "per mil":
                var_data = var_data *1000 # convert to permil
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is isotopic ratio")
   
    # Evaporation            
    elif varname == "evap":
        var_data = Dataset["evap"]
        if units is not None:
            if units == "mm/month":
                var_data = var_data *60*60*24*30 *-1  #mm/month (positive values)
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kg/m²s")
   
    # relative humidity near surface            
    elif varname == "relhum":
        if lev is not None:
            var_data = Dataset["relhum"].sel(lev=lev)
        else:
            var_data = Dataset["relhum"][:,30,:,:]
        if units is not None:
            if units == "%":
                var_data = var_data *100 #%
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is ratio")
                
    # elevation from GEOSP
    elif varname == "elev":
        var_data = Dataset["geosp"] / 9.8  # convert to m**2/s**2 --> m
        if units is not None:
            if units == "km":
                var_data = var_data /1000  #km
                var_data.attrs["units"] = units
            else:
                var_data.attrs["units"] = units
                print("The default units of elevation is used [m]")
                
    # sea land mask
    elif varname == "slm":
        var_data = Dataset["slm"]
        
    # Surface Wind velocity zonal
    elif varname == "u10":
        var_data = Dataset["u10"]
        
    # Surface wind velocity meridoinal
    elif varname == "v10":
        var_data = Dataset["v10"]
        
    elif varname == "wind10":
        var_data = Dataset["wind10"]
    
    # geopotential height at pressure levels
    elif varname == "geopoth": #units in m
        
        var_data = Dataset["geopoth"]
        
        if lev_units is not None:
            
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
        
        if lev is not None:
            
            var_data = var_data.sel(lev=lev)
        
    
    # mean sea level pressure or surface pressure  
    elif varname == "slp":
        if hasattr(Dataset, "aps"):
            var_data = Dataset["aps"]
        elif hasattr(Dataset, "psl"):
            var_data = Dataset["psl"]
        else:    
            var_data = Dataset["slp"]
        
        if units is not None:
            if units == "hPa":
                var_data /= 100
                var_data.attrs["units"] = units
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
    # Vertical verlocity    
    elif varname =="omega":  # Pa/s
        var_data = Dataset[varname]
        
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
            
        if units is not None:
            #Other data must be temperature on pressure levels with the same shape as omega
            if units == "m/s":
                rgas=287.058  # m²/s²K
                g=9.80665
                pa = Dataset["lev"] # must be the same size and shape as omega and t
                rho = pa/(rgas*other_data)
                var_data /= (rho*g)  #Pa/s -->m/s
            elif units == "Pa/s":
                var_data = var_data
                
            else:
                raise ValueError("You have defined incorrect units or its not implemented")
                
        if lev is not None:
            var_data = var_data.sel(lev=lev)
        
        
    # meridional wind at vertical levels    
    elif varname == "v":
        var_data = Dataset[varname]
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
            
        if lev is not None:
            var_data = var_data.sel(lev=lev)
            
    # zonal wind at vertical levels    
    elif varname == "u":
        var_data = Dataset[varname]
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
            
        if lev is not None:
            var_data = var_data.sel(lev=lev)
            
    elif varname == "e/p":
        var_prec = Dataset["aprl"] + Dataset["aprc"]
        var_evap = Dataset["evap"] 
        var_data = var_evap/var_prec
        
    elif varname == "E-P":
        var_prec = Dataset["aprl"] + Dataset["aprc"]
        var_evap = Dataset["evap"] 
        
        if units is not None:
            if units == "mm/month":
                var_prec = var_prec *60*60*24*30  #mm/month (positive values)
                var_evap = var_evap *60*60*24*30*-1  #mm/month (positive values)
                
            else:
                print("Define unit well or the default is kg/m²s")
            
            var_data = var_evap - var_prec
            var_data.attrs["units"] = units
        else:
            var_data = var_evap - var_prec
        
            
    # temperature at vertical levels    
    elif varname == "st":
        var_data = Dataset[varname]
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
            
                 
        if lev is not None:
            var_data = var_data.sel(lev=lev)
            
    # total cloud cover at vertical levels 
    elif varname == "aclcac":
        var_data = Dataset[varname]
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
                 
        if lev is not None:
            var_data = var_data.sel(lev=lev)
            
    # specific humidity at vertical levels         
    elif varname == "q":
        var_data = Dataset[varname]
        
        if units is not None:
            if units =="g/kg":
                var_data = var_data * 1000 #--->g/kg
                
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
                 
        if lev is not None:
            var_data = var_data.sel(lev=lev)
            
    # relativehumidity at vertical levels (pressure)      
    elif varname == "rh":
        var_data = Dataset["relhum"]
        
        if units is not None:
            if units =="%":
                var_data = var_data * 100 #---%
                
        if lev_units is not None:
            var_data = vert_coord_convertion(data=var_data, units=lev_units)
                 
        if lev is not None:
            var_data = var_data.sel(lev=lev)
            
            
    elif varname == "latent heat":
        var_data = Dataset["ahfl"] #positive values
        
    elif varname == "sensible heat":
        var_data = Dataset["ahfs"]
        
    elif varname =="qvi":
        var_data = Dataset["qvi"] # vertical integrated water vapor
        
    elif varname =="apmeb":
        var_data = Dataset["apmeb"] # vertical integrated water vapor
    
    elif varname == "vegetation_frac":
        var_data = Dataset["VGRATCLIM"]
        
    elif varname == "TOA":
        var_data = Dataset["srafl"]
        
    elif varname == "SW_land":
        var_data = Dataset["sofllac"]
        
        
    else:
        print("Check the Varname or it has not been implemeted")
            
    return var_data