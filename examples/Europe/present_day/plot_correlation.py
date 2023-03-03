# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:05:13 2023

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
from pyClimat.plots import plot_correlation

# read the save data


apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.PlateCarree()
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols=3, 
                                                 figsize=(28, 22), subplot_kw={"projection": projection})

plot_correlation(variable="Spearman Coefficients", data=d18op_nao_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                 plot_pvalues=True, pvalue_data=d18op_nao_sig, bottom_labels=False, title="[A] NAO-$\delta^{18}$Op")

plot_correlation(variable="Spearman Coefficients", data=d18op_ea_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax4, fig=fig, plot_pvalues=True, pvalue_data=d18op_ea_sig,
                 title="[D] EA-$\delta^{18}$Op")

plot_correlation(variable="Spearman Coefficients", data=t2m_nao_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax2, fig=fig, plot_pvalues=True, pvalue_data=t2m_nao_sig, bottom_labels=False,
                 left_labels=False, title="[B] NAO-t2m")

plot_correlation(variable="Spearman Coefficients", data=t2m_ea_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax5, fig=fig, plot_pvalues=True, pvalue_data=t2m_ea_sig, left_labels=False,
                 title="[E] EA-t2m")

plot_correlation(variable="Spearman Coefficients", data=prec_nao_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax3, fig=fig, plot_pvalues=True, pvalue_data=prec_nao_sig, bottom_labels=False,
                 left_labels=False, title="[C] NAO-prec")

plot_correlation(variable="Spearman Coefficients", data=prec_ea_sval, units="-", vmax=1,
                 vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax6, fig=fig, plot_pvalues=True, pvalue_data=prec_ea_sig,
                 left_labels=False, title="[F] EA-prec")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.001)
plt.savefig(os.path.join(path_to_plots, "corr_nao_ea_d18op.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.show()