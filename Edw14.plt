Scenarios: # Specify netCDF files and scenario labels
DATA/Edw14M3.nc
"v3"

# Specify cut-off for minor/major reactions
# and net fluxes for inorganic Ox/NOx cycles
Settings:
MCM:     v3.2
cut-off: 0.05, 1.0
cycles:  reduce

Plotting: # Specify plots for pdf

v3: specs/ppb
O3
NO NO2
PAN CH3CHO

v3: specs/ppm
CH4

v3: specs/ppt
HONO
CLNO2

v3: specs
OH
HO2
EMISS

v3: rates/mlc
EMISS-->NO
NO2-->O+NO
O3-->DUMMY
CLNO2-->DUMMY
