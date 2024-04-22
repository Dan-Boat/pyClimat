# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:52:17 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
import cartopy.crs as ccrs
from cycler import cycler

import geocat.comp as gc


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff, extract_transect
from pyClimat.variables import extract_var


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


# define paths
pi_path = "D:/Datasets/Model_output_pst/PI/MONTHLY_MEANS/"
miocene_path= "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"
pliocene_path= "D:/Datasets/Model_output_pst/PLIO/MONTHLY_MEANS/"
holocene_path = "D:/Datasets/Model_output_pst/MH/MONTHLY_MEANS/"
miocene_path = "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"
miocene_ea_high_path = "D:/Datasets/Model_output_pst/a025_dkrz-levante_e5w2.3_t159_MIO_EA_high_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"

echam_dynamics_filename = "1003_1017_1m_mlterm_dynamics.nc"
#echam_filename = "1003_1017_1m_mlterm.nc"

# get all the variables

def mass_weight(pressure_levs):

    g = 9.80665 # gravitational acceleration (m s-2)
    layer_mass_weighting = pressure_levs / g # Layer Mass Weighting
    
    
    layer_mass_weighting = layer_mass_weighting.transpose("lev", "lon", "lat")
    layer_mass_weighting = layer_mass_weighting.rename({"lev":"plev"})
    
    return layer_mass_weighting

def convert_omega(dataset):
    #convert omega to m/s
    rgas=287.058  # m²/s²K
    g=9.80665
    pa = dataset["plev"] # must be the same size and shape as omega and t
    rho = pa/(rgas*dataset["st"])
    omega = dataset.omega/(rho*g)  #Pa/s -->m/s
    
    return omega

def mass_integration(dataset, layer_mass_weighting, omega):
    

    mass_weighted_vapor = dataset.q * layer_mass_weighting # mass weighted 'q'
    mass_weighted_u = dataset.u * layer_mass_weighting # mass weighted 
    mass_weighted_v = dataset.v * layer_mass_weighting # mass weighted 'v'
    mass_weighted_omega = omega * layer_mass_weighting # mass weighted 'omega'
    
    iq = mass_weighted_vapor.sum(dim="plev") # Vertically Integrated Vapor
    iv = mass_weighted_v.sum(dim="plev")
    iu = mass_weighted_u.sum(dim="plev")
    iw = mass_weighted_omega.sum(dim="plev")
    
    du_dx, du_dy = gc.gradient(iu) # (s-1)
    dv_dx, dv_dy = gc.gradient(iv)
    dq_dx, dq_dy = gc.gradient(iq) # (kg m-3)
    
    return iq, iv, iu, iw, du_dx, du_dy, dv_dx, dv_dy, dq_dx, dq_dy

def perform_moisture_diagnostic(control_path, exp_path, max_lon, max_lat, min_lon, 
                                min_lat, time="month", month="JJAS", 
                                filename="1003_1017_1m_mlterm_dynamics.nc"):
    
    control_data = read_from_path(path=control_path, filename=filename)
    
    exp_data = read_from_path(path=exp_path , filename=filename)
    

    control_data_alt = compute_lterm_mean(data=control_data, time=time, month=month)
    exp_data_alt = compute_lterm_mean(data=exp_data, time=time, month=month)
    anomalies_alt = compute_lterm_diff(data_control=control_data, 
                                          data_main=exp_data, time=time, month=month)

        
    region_control = extract_transect(data=control_data_alt, maxlon=max_lon, minlon=min_lon,
                                   maxlat=max_lat, minlat=min_lat)

    region_exp = extract_transect(data=exp_data_alt, maxlon=max_lon, minlon=min_lon,
                                   maxlat=max_lat, minlat=min_lat)

    region_anomalies = extract_transect(data=anomalies_alt, maxlon=max_lon, minlon=min_lon,
                                   maxlat=max_lat, minlat=min_lat)


    exp_delta_p = gc.meteorology.delta_pressure(pressure_lev=region_exp.plev, 
                                                          surface_pressure=region_exp.aps)

    layer_mass_weighting_exp = mass_weight(exp_delta_p)
    
    #omega_pi = convert_omega(extract_wam_pi)
    omega_exp = convert_omega(region_exp)
    omega_anomalies = convert_omega(region_anomalies)

    

    iq, iv, iu, iw, du_dx, du_dy, dv_dx, dv_dy, dq_dx, dq_dy = mass_integration(region_exp, layer_mass_weighting_exp, omega_exp)


    iq_d, iv_d, iu_d, iw_d, du_dx_d, du_dy_d, dv_dx_d, dv_dy_d, dq_dx_d, dq_dy_d = mass_integration(region_anomalies, layer_mass_weighting_exp, omega_anomalies)

    thermodynamic = (-1/100)* 2592000 * (iw *( iq_d * np.add(dq_dx_d.data, dq_dy_d.data)))

    dynamic = (-1)* 2592000 *(iw_d *( iq * np.add(dq_dx.data, dq_dy.data)))

    non_linear = (-1) * 2592000* (iw_d *( iq_d * np.add(dq_dx_d.data, dq_dy_d.data)))
    
    horizontal_advection = 2592000 * (((-iu*dq_dx.data) - (iv*dq_dy.data))/(1e3*9.8))
    
    #vertical_advection = -(iw * (iq * np.add(dq_dx.data, dq_dy.data)))
    
    convergence = 2592000 * (-iq * np.add(du_dx.data, dv_dy.data))/(1e3*9.8)
    
    mfc = horizontal_advection - convergence
    
    
    # thermodynamic = iq_d *( np.add(du_dx.data, dv_dy.data))
    # dynamic = iq * np.add(du_dx_d.data, dv_dy_d.data)
    # non_linear = iq_d * np.add(du_dx_d.data, dv_dy_d.data)
    
    
    p_e_d = (region_anomalies["precip"]-(-1*region_anomalies["evap"]))* 2592000
    p_e = (region_exp["precip"]-(-1*region_exp["evap"])) * 2592000
    
    residual = p_e - horizontal_advection - convergence
    
    # return {"thermodynmic": thermodynamic.sortby(thermodynamic.lon),
    #         "dynamic": dynamic.sortby(dynamic.lon),
    #         "non_linear": non_linear.sortby(non_linear.lon),
    #         "horizontal_advection": horizontal_advection.sortby(horizontal_advection.lon),
    #         "p_e_d": p_e_d.sortby(p_e_d.lon),
    #         "p_e": p_e.sortby(p_e.lon),
    #         "residual": residual.sortby(residual.lon), 
    #         "convergence": convergence.sortby(convergence.lon), 
    #         "mfc": mfc.sortby(mfc.lon)}



    return {"thermodynmic": thermodynamic.sortby(thermodynamic.lon).mean(),
            "dynamic": dynamic.sortby(dynamic.lon).mean(),
            "non_linear": non_linear.sortby(non_linear.lon).mean(),
            "horizontal_advection": horizontal_advection.sortby(horizontal_advection.lon).mean(),
            "p_e_d": p_e_d.sortby(p_e_d.lon).mean(),
            "p_e": p_e.sortby(p_e.lon).mean(),
            "residual": residual.sortby(residual.lon).mean(), 
            "convergence": convergence.sortby(convergence.lon).mean(), 
            "mfc": mfc.sortby(mfc.lon).mean()}
    



moisture_terms_mh = perform_moisture_diagnostic(control_path=pi_path, 
                                             exp_path=holocene_path, max_lon=20,
                                             max_lat=20, min_lon=-10, min_lat=10,
                                             time="month", month="JJAS")


moisture_terms_plio = perform_moisture_diagnostic(control_path=pi_path, 
                                              exp_path=pliocene_path, max_lon=20,
                                              max_lat=20, min_lon=-10, min_lat=10,
                                              time="month", month="JJAS")


# moisture_terms_mio = perform_moisture_diagnostic(control_path=pi_path, 
#                                               exp_path=miocene_path, max_lon=20,
#                                               max_lat=18, min_lon=-10, min_lat=12,
#                                               time="month", month="JJAS")

# moisture_terms_mio_ea_high = perform_moisture_diagnostic(control_path=pi_path, 
#                                               exp_path=miocene_ea_high_path, max_lon=20,
#                                               max_lat=18, min_lon=-10, min_lat=12,
#                                               time="month", month="JJAS")

dTh = "#212F3D"
dD = "#52BE80"
pe = "#CB4335"
nl = "#D4AC0D"
mfc = "magenta"

terms_colors = {"dTH":dTh, "dD":dD, "d(P-E)":pe, "d(NL)":nl}

index = ["MH", "MP"]
columns = ["dTH", "dD", "d(P-E)", "d(NL)"]
df = pd.DataFrame(index=index, columns=columns)

data = [moisture_terms_mh, moisture_terms_plio]

    
for i,x in enumerate(data):
    df["dTH"].loc[index[i]] = float(x.get("thermodynmic").data)
    df["dD"].loc[index[i]] = float(x.get("dynamic").data)
    df["d(P-E)"].loc[index[i]] = float(x.get("p_e_d").data)
    df["d(NL)"].loc[index[i]] = float(x.get("non_linear").data)
    #df["mfc"].loc[index[i]] = float(x.get("mfc").data)
    
    

# plotting

apply_style2(fontsize=24, style="seaborn-talk", linewidth=2.5, usetex=True) 

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 12))

colors = [terms_colors[c] for c in columns] 
mpl.rcParams["axes.prop_cycle"] = cycler("color", colors)

df.plot(kind="bar", legend=True, fontsize=24, width=0.7, ax=ax, rot=45)
ax.set_ylabel("mm/month", fontweight="bold", fontsize=25)
ax.axhline(y=0, linestyle="--", color="black", linewidth=2)
ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1), borderaxespad=0., frameon=True, fontsize=24)
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_fig4.pdf"), format= "pdf", bbox_inches="tight", dpi=600)





moisture_terms_mh = perform_moisture_diagnostic(control_path=pi_path, 
                                             exp_path=holocene_path, max_lon=50,
                                             max_lat=10, min_lon=30, min_lat=-5,
                                             time="month", month="MAM")


moisture_terms_plio = perform_moisture_diagnostic(control_path=pi_path, 
                                              exp_path=pliocene_path, max_lon=50,
                                              max_lat=10, min_lon=30, min_lat=-5,
                                              time="month", month="MAM")


# moisture_terms_mio = perform_moisture_diagnostic(control_path=pi_path, 
#                                               exp_path=miocene_path, max_lon=50,
#                                               max_lat=10, min_lon=30, min_lat=-5,
#                                               time="month", month="MAM")

# moisture_terms_mio_ea_high = perform_moisture_diagnostic(control_path=pi_path, 
#                                               exp_path=miocene_ea_high_path, max_lon=50,
#                                               max_lat=10, min_lon=30, min_lat=-5,
#                                               time="month", month="MAM")


index = ["MH", "MP"]
columns = ["dTH", "dD", "d(P-E)", "d(NL)",]
df = pd.DataFrame(index=index, columns=columns)


data = [moisture_terms_mh, moisture_terms_plio]

    
for i,x in enumerate(data):
    df["dTH"].loc[index[i]] = float(x.get("thermodynmic").data)
    df["dD"].loc[index[i]] = float(x.get("dynamic").data)
    df["d(P-E)"].loc[index[i]] = float(x.get("p_e_d").data)
    df["d(NL)"].loc[index[i]] = float(x.get("non_linear").data)
    #df["mfc"].loc[index[i]] = float(x.get("mfc").data)
    
    

# plotting
apply_style2(fontsize=24, style="seaborn-talk", linewidth=2.5, usetex=True) 

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 12,))
colors = [terms_colors[c] for c in columns] 
mpl.rcParams["axes.prop_cycle"] = cycler("color", colors)

df.plot(kind="bar", legend=True, fontsize=24, width=0.7, ax=ax, rot=45)
ax.set_ylabel("mm/month", fontweight="bold", fontsize=25)
ax.axhline(y=0, linestyle="--", color="black", linewidth=2)
ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1), borderaxespad=0., frameon=True, fontsize=24)

plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_fig5.pdf"), format= "pdf", bbox_inches="tight", dpi=600)






            
# projection = ccrs.Robinson(central_longitude=0, globe=None)

# pprojection_trans = ccrs.PlateCarree()

# fig,(ax1, ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=4, figsize=(18,13), subplot_kw={"projection":projection})

# plot_annual_mean(variable="P-E", data_alt=moisture_terms_mh.get("p_e"), ax=ax1,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(a) P-E", 
#                 domain="West Africa", cbar_pos= [0.35, 0.01, 0.25, 0.02], center=True,
#                 orientation="horizontal")

# plot_annual_mean(variable="Precipitation", data_alt=moisture_terms_mh.get("horizontal_advection"), ax=ax2,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) H. advection", 
#                 domain="West Africa", center=True)


# plot_annual_mean(variable="Precipitation", data_alt=moisture_terms_mh.get("convergence"), ax=ax3,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) Convergence", 
#                 domain="West Africa", center=True)

# plot_annual_mean(variable="Precipitation", data_alt=moisture_terms_mh.get("residual"), ax=ax4,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(d) residual", 
#                 domain="West Africa", center=True)




# fig,(ax1, ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=4, figsize=(18,13), subplot_kw={"projection":projection})

# plot_annual_mean(variable="P-E", data_alt=moisture_terms_mh.get("p_e_d"), ax=ax1,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(a) P'-E'", 
#                 domain="West Africa", cbar_pos= [0.35, 0.01, 0.25, 0.02], center=True,
#                 orientation="horizontal")

# plot_annual_mean(variable="Precipitation", data_alt=moisture_terms_mh.get("thermodynmic"), ax=ax2,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) Thermodynamic", 
#                 domain="West Africa", center=True)


# plot_annual_mean(variable="Precipitation", data_alt=moisture_terms_mh.get("dynamic"), ax=ax3,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) dynamic", 
#                 domain="West Africa", center=True)

# plot_annual_mean(variable="Precipitation", data_alt=moisture_terms_mh.get("non_linear"), ax=ax4,
#                  cmap="RdBu", units="mm/month", vmax=100, vmin=-100, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(d) NL", 
#                 domain="West Africa", center=True)

# # plotting
# apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
# projection = ccrs.Robinson(central_longitude=0, globe=None)

# pprojection_trans = ccrs.PlateCarree()

# fig,(ax1, ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=4, figsize=(30,13), subplot_kw={"projection":projection})


# (advection*60*60*24*30).sortby(advection.lon).plot(ax=ax1)
# (convergence*60*60*24*30).sortby(convergence.lon).plot(ax=ax2)
# (extract_wam["precip"]*60*60*24*30).sortby(extract_wam.lon).plot(ax=ax3)
# (extract_wam["evap"]*60*60*24*30).sortby(extract_wam.lon).plot(ax=ax4)





