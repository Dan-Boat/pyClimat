#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:53:11 2021

@author: dboateng
This routine contains funtional utilities required in the other modules
"""
def vert_coord_convertion(data, units):
    """
    

    Parameters
    ----------
    data : TYPE: Dataarray (MxNxV)
        DESCRIPTION. Dataarray on vertical coordinates with Pa as default unit
    units : TYPE: str
        DESCRIPTION. The required unit for conversion eg. hPa

    Raises
    ------
    ValueError
        DESCRIPTION. When other units is defined aside hPa 

    Returns
    -------
    data : TYPE: Datarray 
        DESCRIPTION. Data vertical units converted

    """
    if hasattr(data, "mlev"):
        data = data.rename({"mlev":"lev"})
    elif hasattr(data, "plev"):
        data = data.rename({"plev":"lev"})
        
    if units =="hPa":
        data["lev"] = data.lev/100  #Pa -->hPa
        data["lev"].attrs["units"] = "hPa"
    else:
        raise ValueError("The units may be incorrect or not implemented yet!")
    
    return data


