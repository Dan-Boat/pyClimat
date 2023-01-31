# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:35:52 2023

@author: dboateng

This module contains all the statistic functions used in pyClimat (eg. EOF, lingress, t-test)
"""

# To do list
#1. Implement EOF class that has the various methods
#2. Implement statistical testing class for xarray datasets


# importing packages 
import xarray as xr
import os
import pandas as pd
import numpy as np
from scipy import stats
from eofs.xarray import Eof
import scipy as sc
from sklearn.decomposition import PCA

def _get_month(npdatetime64):
    
     """
     Returns the month for a given npdatetime64 object, 1 for January, 2 for
     February, ...
     """
     month =  npdatetime64.astype('datetime64[M]').astype(int) % 12 + 1
     
     return month


class EOF_standard():
    
    def __init__(self, data, weights=True, standardize=True,
                 extract_region=True, extract_season=True, neofs=4,):
        
        self.data = data
        
        if not isinstance(self.data, xr.DataArray):
            raise TypeError("The X field must be a datarray object or type")
            
        self.weights = weights
        self.standardize = standardize
        self.extract_region = extract_region
        self.extract_season = extract_season
        self.neofs = neofs
    
    
    
    
    def select_time_and_region(self, maxlon, minlon, maxlat, minlat, time="season", 
                               month=None, season=None):
        
        if self.extract_region:
            if hasattr(self.data, "longitude"):
                self.data = self.data.rename({"longitude":"lon", "latitude":"lat"})
                
            
            data = self.data.assign_coords({"lon": (((self.data.lon + 180) % 360) - 180)})
    
            data = data.where((data.lat >= minlat) & (data.lat <= maxlat), drop=True)
            data = data.where((data.lon >= minlon) & (data.lon <= maxlon), drop=True)
            
            self.data = data
            
            
        if self.extract_season:
            
            if time =="season":
                data_group = self.data.groupby("time.season")
                
                if season is not None:
                    data = data_group[season]
                    
                else:
                    data = data_group["DJF"]  # set to NH Winter as default
                    
            elif time == "month":
                
                data_group = data.groupby("time.month")
                
                if isinstance(month, int):
                    data = data_group.get(month)
                
                elif month == "ONDJFM":
                    data = data.sel(time=data.time.dt.month.isin([1, 2, 3, 10, 11, 12]))
                    
                elif month == "AMJJAS":
                    data = data.sel(time=data.time.dt.month.isin([4, 5, 6, 7, 8, 9]))
                    
                else:
                    raise ValueError("The define month parameter is not recognized")
        
        self.data = data                
        
    
    def calculate_anomalies(self, if_return=False):
        # remove monthly cyle
        monthly_means = self.data.groupby("time.month").mean(dim="time")
        
        
        group = self.data.groupby("time.month")
        
        
        anomalies = group.apply(
            lambda x: x - monthly_means.sel(month=_get_month(x[0].time.values))
        )
        #self.anomalies = self.data - monthly_means
        
        self.anomalies = anomalies.drop("month")
        
        if self.standardize:
            self.anomalies /= self.anomalies.std(dim="time")
            
        if self.weights:
            self.anomalies = self.anomalies * np.sqrt(np.abs(np.cos(self.anomalies.lat*np.pi/180))) 
            
        
        if True in self.anomalies.isnull():
            self.anomalies = self.anomalies.dropna(dim="time")
            
            
        if if_return:
            
            return self.anomalies
        
        
    
    
    
    def eof_solver(self, method="Eof", apply_varimax=False, apply_promax=False,):
        """
        *****************************************
        n_rot : int
        Number of modes to be rotated.
            power : int
                Defines the power of Promax rotation. Choosing ``power=1`` equals
                a Varimax solution (the default is 1).
            max_iter : int
                Number of maximal iterations for obtaining the rotation matrix
                (the default is 1000).
            rtol : float
                Relative tolerance to be achieved for early stopping the iteration
                process (the default is 1e-8).
        *****************************************
        Parameters
        ----------
        method : TYPE, optional
            DESCRIPTION. The default is "Eof".
        apply_varimax : TYPE, optional
            DESCRIPTION. The default is False.
        apply_promax : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        
        self.method = method
        self.apply_varimax = apply_varimax
        self.apply_promax = apply_promax
        
        
        if self.method == "Eof":
            self.solver = Eof(self.anomalies)   
            
        elif self.method =="xeofs":
            
            from xeofs.xarray import EOF
            
            self.model = EOF(self.anomalies, n_modes=self.neofs, dim="time", norm=False)
            
            self.model.solve()
            
            if self.apply_varimax:
                from xeofs.xarray import Rotator
                
                self.rot_varimax = Rotator(self.model, n_rot=50, power=1)
                
            if self.apply_promax:
                
                from xeofs.xarray import Rotator
                
                self.rot_promax = Rotator(self.model, n_rot=50, power=4)
                
            
            
    def eofs(self, eofscaling=0):
        """
        *eofscaling for Eof*
            Sets the scaling of the EOFs. The following values are
            accepted:

            * *0* : Un-scaled EOFs (default).
            * *1* : EOFs are divided by the square-root of their
              eigenvalues.
            * *2* : EOFs are multiplied by the square-root of their
              eigenvalues.
              
         *eofscaling : [0, 1, 2] for xeofs
         
         EOFs are scaled (i) to be orthonormal (``scaling=0``), (ii) by the
         square root of the eigenvalues (``scaling=1``) or (iii) by the
         singular values (``scaling=2``). In case no weights were applied,
         scaling by the singular values results in the EOFs having the
         unit of the input data (the default is 0).

        *neofs*
            Number of EOFs to return. Defaults to all EOFs. If the
            number of EOFs requested is more than the number that are
            available, then all available EOFs will be returned.

        Parameters
        ----------
        eofscaling : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        """
        
        if self.method == "Eof":
            self.eofs = self.solver.eofs(eofscaling=eofscaling, neofs=self.neofs)
            
                
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.eofs = self.rot_varimax.eofs(scaling=eofscaling)   
                
            elif self.apply_promax:
                self.eofs = self.rot_promax.eofs(scaling=eofscaling) 
            
            
            else:
                self.eofs = self.model.eofs(scaling=eofscaling)
                
        self.eofs = self.eofs.sortby(self.eofs.lon)    # fix the change of +-180 changes
        
        return self.eofs
    
    
    def pcs(self, pscaling=0):
        
        if self.method == "Eof":
            self.pcs = self.solver.pcs(pcscaling=pscaling, npcs=self.neofs)
                
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.pcs = self.rot_varimax.pcs(scaling=pscaling)   
                
            elif self.apply_promax:
                self.pcs = self.rot_promax.pcs(scaling=pscaling) 
            
            else:
                self.pcs = self.model.pcs(scaling=pscaling)
         
                
        return self.pcs.to_pandas()
        
    
    def eofs_as_correlation(self):
        
        if self.method == "Eof":
            self.eofsAsCorrelation = self.solver.eofsAsCorrelation(neofs=self.neofs)
            
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.eofsAsCorrelation = self.rot_varimax.eofs_as_correlation()  
                
            elif self.apply_promax:
                self.eofsAsCorrelation = self.rot_promax.eofs_as_correlation() 
            
            else:
                self.eofsAsCorrelation = self.model.eofs_as_correlation()
        
        
        self.eofsAsCorrelation = self.eofsAsCorrelation.sortby(self.eofsAsCorrelation.lon)
        
        return self.eofsAsCorrelation
    
    
    
    def explained_variance(self):
        
        if self.method == "Eof":
            
            self.explained_variance = self.solver.totalAnomalyVariance()
        
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                
                self.explained_variance = self.rot_varimax.explained_variance()
                
            elif self.apply_promax:
                
                self.explained_variance = self.rot_promax.explained_variance()
                
            else: 
                self.explained_variance = self.model.explained_variance()
        
        return self.explained_variance
    
    
    
    def explained_variance_ratio(self):
        
        if self.method == "Eof":
            
            self.explained_variance_ratio = self.solver.varianceFraction(neigs=self.neofs)
        
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                
                self.explained_variance_ratio = self.rot_varimax.explained_variance_ratio()
                
            elif self.apply_promax:
                
                self.explained_variance_ratio = self.rot_promax.explained_variance_ratio()
                
            else: 
                self.explained_variance_ratio = self.model.explained_variance_ratio()
        
        return self.explained_variance_ratio
    
    def reconstuct_X(self, mode=None):
        
        if self.method == "Eof":
            X_field = self.solver.reconstructedField(neofs=self.neofs)
            
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                
                if mode is not None:
                    X_field = self.rot_varimax.reconstruct_X(mode=mode)
                    
                else:
                    X_field = self.rot_varimax.reconstruct_X(mode=1)
                    
            elif self.apply_promax:
                if mode is not None:
                    X_field = self.rot_promax.reconstruct_X(mode=mode)
                    
                else:
                    X_field = self.rot_promax.reconstruct_X(mode=1)
        
        pass
    
    def project_X_onto_eofs(self, data, eofscaling=0, weighted=False):
        
        if self.method == "Eof":
            self.pcs_projected = self.solver.projectField(array=data, neofs=self.noefs, 
                                                          eofscaling=eofscaling)
        
        
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.pcs_projected = self.rot_varimax
                
                
        pass
    


    