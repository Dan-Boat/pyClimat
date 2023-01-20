#!/bin/bash
### This script calculates the long-term means of the downloaded PMIP4 files and also merge time into one file
#
# @ Daniel Boateng 
# The script is based on cdo for the computation and would be analysed in python
# The data were downloaded from the cmip6 host databases (dkrz) and also from the /poo/data/CMIP on levante
#
#Defining paths for the experiments

OPT=4   # 1 for MH, 2 for LGM, 3 for mPLIO, and 4 for PI
mPATH="/esd/esd02/data/climate_models/PMIP4"
MH_FOLDER="MidHolocene"
LGM_FOLDER="lgm"
PLIO_FOLDER="MidPliocene"
PI_FOLDER="Pre-Industrial"


MH_MODELS=("AWI.AWI-ESM-1-1-LR"
"CESM2"
"EC-Earth3-LR"
"IPSL-CM6A-LR"
"MIROC-ES2L"
"MPI-ESM1-2-LR"
"GISS-E2-1-G"
"HadGEM3-GC31-LL")


LGM_MODELS=("AWI-ESM-1-1-LR"
"CESM2-WACCM-FV2"
"INM-CM4-8"
"MIROC-ES2L"
"MPI-ESMI-2-LR"
)

PLIO_MODELS=("CESM2"
"EC-Earth3-LR"
"GISS-E2-1-G"
"HadGEM3-GC31-LL"
"IPSL-CM6A-LR"
"NorESM1-F"
)

PI_MODELS=("AWI-ESM-1-1-LR"
"CESM2"
"EC-Earth3-LR"
"IPSL-CM6A-LR"
"MIROC-ES2L"
"MPI-ESM1-2-LR"
"GISS-E2-1-G"
"HadGEM3-GC31-LL"
"CESM2-WACCM-FV2"
"INM-CM4-8"
"NorESM1-F")


module load cdo_1.9.9
# loog throught the models and extract all the nc files, then merge then store there, take the merge file and calculate the lterm mean
if ((${OPT} == "1")); then
# long term for MH
# =========================================================================
for i in {0..7}; do
    
    infile_pr=${mPATH}/${MH_FOLDER}/${MH_MODELS[$i]}/pr/pr_*.nc

    cdo mergetime ${infile_pr} ${mPATH}/${MH_FOLDER}/postprocessed/${MH_MODELS[$i]}/pr_monthly.nc
    cdo ymonmean ${mPATH}/${MH_FOLDER}/postprocessed/${MH_MODELS[$i]}/pr_monthly.nc ${mPATH}/${MH_FOLDER}/postprocessed/${MH_MODELS[$i]}/pr_1m_lterm.nc

    infile_tas=${mPATH}/${MH_FOLDER}/${MH_MODELS[$i]}/tas/tas_*.nc

    cdo mergetime ${infile_tas} ${mPATH}/${MH_FOLDER}/postprocessed/${MH_MODELS[$i]}/tas_monthly.nc
    cdo ymonmean ${mPATH}/${MH_FOLDER}/postprocessed/${MH_MODELS[$i]}/tas_monthly.nc ${mPATH}/${MH_FOLDER}/postprocessed/${MH_MODELS[$i]}/tas_1m_lterm.nc
    
done
elif ((${OPT} == "2")); then
# long term for LGM
# =========================================================================

for i in {0..5}; do
    
    infile_pr=${mPATH}/${LGM_FOLDER}/${LGM_MODELS[$i]}/pr/pr_*.nc

    cdo mergetime ${infile_pr} ${mPATH}/${LGM_FOLDER}/postprocessed/${LGM_MODELS[$i]}/pr_monthly.nc
    cdo ymonmean ${mPATH}/${LGM_FOLDER}/postprocessed/${LGM_MODELS[$i]}/pr_monthly.nc ${mPATH}/${LGM_FOLDER}/postprocessed/${LGM_MODELS[$i]}/pr_1m_lterm.nc

    infile_tas=${mPATH}/${LGM_FOLDER}/${LGM_MODELS[$i]}/tas/tas_*.nc

    cdo mergetime ${infile_tas} ${mPATH}/${LGM_FOLDER}/postprocessed/${LGM_MODELS[$i]}/tas_monthly.nc
    cdo ymonmean ${mPATH}/${LGM_FOLDER}/postprocessed/${LGM_MODELS[$i]}/tas_monthly.nc ${mPATH}/${LGM_FOLDER}/postprocessed/${LGM_MODELS[$i]}/tas_1m_lterm.nc
    
done
elif ((${OPT} == "3")); then
# long term fpr PLIO
# =========================================================================
for i in {0..5}; do
    
    infile_pr=${mPATH}/${PLIO_FOLDER}/${PLIO_MODELS[$i]}/pr/pr_*.nc

    cdo mergetime ${infile_pr} ${mPATH}/${PLIO_FOLDER}/postprocessed/${PLIO_MODELS[$i]}/pr_monthly.nc
    cdo ymonmean ${mPATH}/${PLIO_FOLDER}/postprocessed/${PLIO_MODELS[$i]}/pr_monthly.nc ${mPATH}/${PLIO_FOLDER}/postprocessed/${PLIO_MODELS[$i]}/pr_1m_lterm.nc

    infile_tas=${mPATH}/${PLIO_FOLDER}/${PLIO_MODELS[$i]}/tas/tas_*.nc

    cdo mergetime ${infile_tas} ${mPATH}/${PLIO_FOLDER}/postprocessed/${PLIO_MODELS[$i]}/tas_monthly.nc
    cdo ymonmean ${mPATH}/${PLIO_FOLDER}/postprocessed/${PLIO_MODELS[$i]}/tas_monthly.nc ${mPATH}/${PLIO_FOLDER}/postprocessed/${PLIO_MODELS[$i]}/tas_1m_lterm.nc
    
done
elif ((${OPT} == "4")); then
# long term for PI
# =========================================================================
for i in {0..11}; do
    
    infile_pr=${mPATH}/${PI_FOLDER}/${PI_MODELS[$i]}/pr/pr_*.nc

    cdo mergetime ${infile_pr} ${mPATH}/${PI_FOLDER}/postprocessed/${PI_MODELS[$i]}/pr_monthly.nc
    cdo ymonmean ${mPATH}/${PI_FOLDER}/postprocessed/${PI_MODELS[$i]}/pr_monthly.nc ${mPATH}/${PI_FOLDER}/postprocessed/${PI_MODELS[$i]}/pr_1m_lterm.nc

    infile_tas=${mPATH}/${PI_FOLDER}/${PI_MODELS[$i]}/tas/tas_*.nc

    cdo mergetime ${infile_tas} ${mPATH}/${PI_FOLDER}/postprocessed/${PI_MODELS[$i]}/tas_monthly.nc
    cdo ymonmean ${mPATH}/${PI_FOLDER}/postprocessed/${PI_MODELS[$i]}/tas_monthly.nc ${mPATH}/${PI_FOLDER}/postprocessed/${PI_MODELS[$i]}/tas_1m_lterm.nc
    
done
fi