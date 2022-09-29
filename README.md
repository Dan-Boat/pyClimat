
<h1 align="center">Welcome to REAME pyClimat 👋</h1>

<p align="center">
<a href="https://pypi.python.org/pypi/py" target="_blank">
  <img src="https://img.shields.io/pypi/v/pyClimat.svg" alt="PyPi">
</a>
<a href="https://pypi.org/project/pyClimat" target="_blank">
  <img src="https://img.shields.io/pypi/pyversions/pyClimat" alt="PyPI - Python Version">
</a>
</h1>

Climat python package for analysising GCM model output and visualization. The package is written in a function based 
(would be pivoted to OPP style in future development). The analysis module features climate variable extraction 
and estimation of statistical long-term means and difference. Statistical tools like PCA or EOF analysis are included for specific 
estimates and many other classical methods like testing, OLS estiates, etc. 

## installation 

The easy way is to use `pip intall pyClimat` (but would require some dependencies)
The following packages must be installed: cartopy and xarray. If you failed to compile this in your environment,
kindy raise an issue on that so I rebuid the dist for all systems. The stable verison should work on UNIX platforms for now

Alternatively, the package can be installed with -e flag in edit mode

## Documentation 

The docs folder [docs](./docs/) contains all the dist files for the documentation compilation: 
 `cd docs | make html`
The direct stable link would be shared shortly

## Examples

This package was adopted for all the visualization in the research study by Boateng et. al 2022 ( Impacts of surface uplift on regional climate)
The scripts can be located in the folder [Alps](./examples/Alps/)


- ⚡ Fun fact **Happy coding and contact me if you have issues with the package**

