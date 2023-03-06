# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:35:52 2023

@author: dboateng

This module contains all the statistical functions used in pyClimat (eg. EOF, lingress, t-test,
                                                                     causality, cross-correlation)
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
from xarray import DataArray
from eofs.xarray import Eof
import scipy as sc
from sklearn.decomposition import PCA
from statsmodels.tsa.stattools import acf, adfuller, ccf, grangercausalitytests
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates 
from sklearn.preprocessing import StandardScaler
import scipy.special as special
from scipy.ndimage import uniform_filter



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
                    
                elif month == "JA": # July and August for the SNAO
                    data = data.sel(time=data.time.dt.month.isin([7, 8]))
                    
                else:
                    raise ValueError("The define month parameter is not recognized")
        
        self.data = data                
        
    
    def calculate_anomalies(self, if_return=False, monthly_anomalies=True):
        
        if monthly_anomalies:
            # remove monthly cyle
            monthly_means = self.data.groupby("time.month").mean(dim="time")
            
            
            group = self.data.groupby("time.month")
            
            
            anomalies = group.apply(
                lambda x: x - monthly_means.sel(month=_get_month(x[0].time.values))
            )
            
            self.anomalies = anomalies.drop("month")
            
            
        else:
            
            self.anomalies = self.data - self.data.mean(dim="time")
        
        
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
            
            self.model = EOF(self.anomalies, n_modes=self.neofs, dim="time", norm=True)
            
            self.model.solve()
            
            if self.apply_varimax:
                from xeofs.xarray import Rotator
                
                self.rot_varimax = Rotator(self.model, n_rot=self.neofs, power=1)
                
            if self.apply_promax:
                
                from xeofs.xarray import Rotator
                
                self.rot_promax = Rotator(self.model, n_rot=self.neofs, power=4)
                
            
            
    def eofs(self):
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
            self.eofs = self.solver.eofsAsCovariance(pcscaling=1, neofs=self.neofs)
            
                
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.eofs = self.rot_varimax.eofs(scaling=2)   
                
            elif self.apply_promax:
                self.eofs = self.rot_promax.eofs(scaling=2) 
            
            else:
                self.eofs = self.model.eofs(scaling=2)
                
        self.eofs = self.eofs.sortby(self.eofs.lon)    # fix the change of +-180 changes
        
        return self.eofs
    
    
    def pcs(self):
        
        if self.method == "Eof":
            self.pcs = self.solver.pcs(pcscaling=1, npcs=self.neofs)
                
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.pcs = self.rot_varimax.pcs(scaling=1)   
                
            elif self.apply_promax:
                self.pcs = self.rot_promax.pcs(scaling=1) 
            
            else:
                self.pcs = self.model.pcs(scaling=1)
                
            self.pcs = self.pcs / self.pcs.std(dim="time")
         
                
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
        
        
        self.corr = self.eofsAsCorrelation[0]
        self.pvals = self.eofsAsCorrelation[1]
        
        self.corr = self.corr.sortby(self.corr.lon)
        self.pvals = self.pvals.sortby(self.pvals.lon)
        
        return self.corr, self.pvals
    
    
    
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
        
        return self.explained_variance.to_pandas()
    
    
    
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
        
        return self.explained_variance_ratio.to_pandas()
    
    def reconstuct_X(self, mode=1):
        
        if self.method == "Eof":
            X_field = self.solver.reconstructedField(neofs=self.neofs)
            
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                
                X_field = self.rot_varimax.reconstruct_X(mode=mode)
                    
                    
            elif self.apply_promax:
                
                X_field = self.rot_promax.reconstruct_X(mode=mode)
    
            else:
                
                X_field = self.model.reconstruct_X(mode=mode)

    
    def project_X_onto_eofs(self, data, eofscaling=0):
        
        if self.method == "Eof":
            self.pcs_projected = self.solver.projectField(array=data, neofs=self.noefs, 
                                                          eofscaling=eofscaling, weighted=False)
        
        elif self.method == "xeofs":
            
            if self.apply_varimax:
                self.pcs_projected = self.rot_varimax.project_onto_eofs(X=data, scaling=eofscaling)
                
            
            elif self.apply_promax:
                
                self.pcs_projected = self.rot_promax.project_onto_eofs(X=data, scaling=eofscaling)
                
            else:
                self.pcs_projected = self.model.project_onto_eofs(X=data, scaling=eofscaling)
            
    

def StackArray(x,dim):
	'''return stacked array with only one dimension left
	INPUTS:
	   x  : xarray.DataArray or Dataset to be stacked
	   dim: sole dimension to remain after stacking
	OUTPUTS:
	   stacked: same as x, but stacked
	'''
	dims = []
	for d in x.dims:
		if d != dim:
			dims.append(d)
	return x.stack(stacked=dims)



def ComputeCorr(i, x, y, method="Spearmanr"):
    
    if x.shape == y.shape:
        
        if method == "Spearmanr":
            sloc,ploc = stats.spearmanr(x.isel(stacked=i), y.isel(stacked=i))
            
        elif method == "pearsonr":
            sloc,ploc = stats.pearsonr(x.isel(stacked=i), y.isel(stacked=i))
            
    else:
        
        if method == "Spearmanr":
            sloc,ploc = stats.spearmanr(x.isel(stacked=i), y.isel(stacked=0))
            
        elif method == "pearsonr":
            sloc,ploc = stats.pearsonr(x.isel(stacked=i), y.isel(stacked=0))
        
    return sloc, ploc
    
def StatCorr(x,y,dim=None, return_sig=True, sig=0.1):
    if len(x.time) != len(y.time):
        x = x.drop_duplicates(dim="time")
        y = y.drop_duplicates(dim="time")
        
        
    if len(y.dims) ==1 or dim is None:
        
        sy = y.expand_dims(stacked=[0])
    else:
        sy = StackArray(x=y, dim=dim)
    
    if dim is None or len(x.dims) == 1:
        sx = x.expand_dims(stacked=[0])
    else:
        sx = StackArray(x,dim)
        
    nspace = len(sx.stacked)
    sval, pval = np.zeros(sx.stacked.shape), np.zeros(sx.stacked.shape)
    for i in range(nspace):
        sval[i], pval[i] = ComputeCorr(i, sx, sy)
        
    if nspace > 1:
        pvalx = DataArray(pval, coords=[sx.stacked],name='pval').unstack('stacked')
        svalx = DataArray(sval,coords=[sx.stacked],name='sval').unstack('stacked')
        svalx, pvalx = svalx.sortby("lon"), pvalx.sortby("lon")
        
    else:
        svalx, pvalx = pval[0], sval[0]
        
    sig_loc  = xr.where(pvalx < sig, pvalx, pvalx*np.nan)
    
    #sig_loc = sig_loc.sortby("lon")
    
    
    
    
    if return_sig:
        return svalx, pvalx, sig_loc
        
    else:
        return svalx, pvalx

def StatCrossCorr(x,y, dim=None, plot=False, sample_rate=3, apply_standardize=False):
    if x.shape != y.shape:
        raise ValueError("Both data should have the same time axis")
        
    sy = y.expand_dims(stacked=[0]).isel(stacked=0)
    sx = x.expand_dims(stacked=[0]).isel(stacked=0)
    
    if apply_standardize == True:
        sx = sx.values.reshape(-1,1)
        sy = sy.values.reshape(-1,1)
        
        scaler = StandardScaler()
        
        scaler_x = scaler.fit(sx)
        scaler_y = scaler.fit(sy)
        
        sx = scaler_x.transform(sx)
        sy = scaler_y.transform(sy)
    
    coefs = ccf(sx, sy, adjusted=False)
    
    # remove padding and reverse order 
    #coefs = coefs[0:(len(x) + 1)][::-1]
    
    if plot:
        plt.subplots(figsize=(15, 4))
        plt.ylabel('Correlation Coefficient')
        plt.xlabel('Lag [time] $h$')
        plt.stem(np.linspace(0, len(x), sample_rate*len(x)), coefs)
        peak = np.argmax(coefs/sample_rate)
        plt.axvline(peak, c='r')
        print('Maximum at:', peak, '(in $s$)')
        plt.show()
        
        
    return coefs
        
class GrangerCausality():
    """
    Note that this class is for spatial data at the monemt (station data would be implemented 
                                                            in a differnt class or adjust this
                                                            one)
    """
    
    def __init__(self, maxlag=10, test="params_ftest"):
        self.maxlag = maxlag
        self.test = test
        
    def prepare_data(self, X, Y, Z=None, dim="time", ):
        
        if isinstance(Y, xr.DataArray):
            if len(Y.dims) ==1 or dim is None:
                sy = Y.expand_dims(stacked=[0])
            else:
                sy = StackArray(x=Y, dim=dim)
        
        if isinstance(X, xr.DataArray):
            
            if dim is None or len(X.dims) == 1:
                sx = X.expand_dims(stacked=[0])
            else:
                sx = StackArray(x=X,dim=dim)
        
        if Z is not None:
            if isinstance(Z, xr.DataArray):
                if dim is None or len(Z.dims) == 1:
                    sz = Z.expand_dims(stacked=[0])
                else:
                    sz = StackArray(x=Z,dim=dim)
                    
                return sx, sy, sz
            else:
                raise ValueError("The alternative variable should be dataaray")
        else:
                
            return sx, sy
    
    def compute_granger(self, i, x, y, z=None, apply_standardize=False, 
                             interchange=False):
        
        if x.shape == y.shape:
            xi = x.isel(stacked=i)
            yi = y.isel(stacked=i)
        elif x.shape > y.shape:
            xi = x.isel(stacked=i)
            yi = y.isel(stacked=0)
        
        else:
            xi = x.isel(stacked=0)
            yi = y.isel(stacked=i)
        
        if z is not None:
            if len(z.stacked) > 1:
                zi = z.isel(stacked=i)
            else:
                zi = z.isel(stacked=0)
                
                
        if apply_standardize == True:
            scaler = StandardScaler()
            
            xi = xi.values.reshape(-1,1)
            scaler_x = scaler.fit(xi)
            xi = scaler_x.transform(xi)
            
            
            yi = yi.values.reshape(-1,1)
            scaler_y = scaler.fit(yi)
            yi = scaler_y.transform(yi)
            
            if z is not None:
                zi = zi.values.reshape(-1,1)
                scaler_z = scaler.fit(zi)
                zi = scaler_z.transform(zi)
                
            
        if interchange:
            if z is not None:
                data = np.column_stack([yi, xi, zi])  # for checking the reverse order of causing
                
            else:
                data = np.column_stack([yi, xi])
        else:
            if z is not None:
                data = np.column_stack([xi, yi, zi])
                
            else:
                data = np.column_stack([xi, yi])
                
                
        # if alternative is provided, then the VAR model should be used directly
        if z is not None:
            
            from statsmodels.tsa.api import VAR
            model = VAR(data)
            model = model.fit(maxlags=self.maxlag, method="bic")
            stats = model.test_causality(caused=0, causing=[1,2], signif=0.05)
        else:
            
            stats = grangercausalitytests(x=data, maxlag=self.maxlag, verbose=False)
        
        return stats   
        
            
            
            
                
    def perform_granger_test(self, X, Y, Z=None, dim="time", apply_standardize=False, 
                             interchange=False):
        
        if Z is not None:
            x,y,z = self.prepare_data(X, Y, Z, dim)
            
        else:
            x,y = self.prepare_data(X, Y, Z, dim)
            z = None
        
        if len(X.shape) != 1:
            nspace = len(x.stacked)
            pval = np.zeros(x.stacked.shape)
        else:
            nspace = len(y.stacked)
            pval = np.zeros(y.stacked.shape)
        
        
        for i in range(nspace):
                
            stats = self.compute_granger(i, x, y, z, apply_standardize, interchange)
            
            if Z is not None:
                pval[i] = stats.pvalue
                
            else:
                pvalues = [round(stats[i+1][0][self.test][1], 4) for i in range(self.maxlag)]
                pvalue = np.min(pvalues)
                pval[i] = pvalue
            
        if len(X.shape) > 1:
            pvalx = DataArray(pval, coords=[x.stacked],name='pval').unstack('stacked')
            pvalx = pvalx.sortby("lon")
            
        elif len(Y.shape) > 1:
            pvalx = DataArray(pval, coords=[y.stacked],name='pval').unstack('stacked')
            pvalx = pvalx.sortby("lon")
        else:
            pvalx = pval[0]
            
        return pvalx        
    
    

def sliding_correlation(a,b,W, method="df_corr", sig=10, plot=False):
    
    
    if method == "df_corr":
        
        df_coef = a.rolling(W).corr(b)
        deg_of_free = W - 2
        t_squared = df_coef**2 * (deg_of_free / ((1.0 - df_coef) * (1.0 + df_coef)))
        pval = special.betainc(0.5*deg_of_free, 0.5, deg_of_free/(deg_of_free + t_squared))
        
        if sig is not None:
            coef_sig = df_coef[pval < sig/100].max()
            
    else:
        # a,b are input arrays; W is window length
        time = a.index
        
        a = a.to_numpy(copy=False)
        b = b.to_numpy(copy=False)
        
        df_coef = pd.DataFrame(index=time, columns= ["coef"])
        am = uniform_filter(a.astype(float),W)
        bm = uniform_filter(b.astype(float),W)
    
        amc = am[W//2:-W//2+1]
        bmc = bm[W//2:-W//2+1]
    
        da = a[:,None]-amc
        db = b[:,None]-bmc
    
        # Get sliding mask of valid windows
        m,n = da.shape
        mask1 = np.arange(m)[:,None] >= np.arange(n)
        mask2 = np.arange(m)[:,None] < np.arange(n)+W
        mask = mask1 & mask2
        dam = (da*mask)
        dbm = (db*mask)
    
        ssAs = np.einsum('ij,ij->j',dam,dam)
        ssBs = np.einsum('ij,ij->j',dbm,dbm)
        D = np.einsum('ij,ij->j',dam,dbm)
        coeff = D/np.sqrt(ssAs*ssBs)
    
        n = W
        ab = n/2 - 1
        pval = 2*special.btdtr(ab, ab, 0.5*(1 - abs(np.float64(coeff))))
        df_coef.iloc[W-1:, 0] = coeff
        
        if sig is not None:
            coef_sig = np.min(coeff[pval < sig/100])
            
    
    if plot:
        fig, ax = plt.subplots(figsize=(15, 4))
        ax.set_ylabel('Correlation Coefficient', fontweight="bold", fontsize=20)
        ax.plot(df_coef, color="black")
        ax.xaxis.set_major_locator(YearLocator(5))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.axhline(y=np.abs(coef_sig), linestyle="--", color="red", linewidth=2)
        ax.axhline(y=0, linestyle="--", color="grey", linewidth=2)
        ax.axhline(y=np.abs(coef_sig)*-1, linestyle="--", color="red", linewidth=2)
        plt.show()
            
    if sig is not None:
        return df_coef, np.abs(coef_sig)
    
    else:
        return df_coef
    