Scenarios: # Specify netCDF files and scenario labels
DATA/Edw14M3.nc
"v3"

# Specify cut-off for minor/major reactions
# and net fluxes for inorganic Ox/NOx cycles
Settings:
cut-off: 0.05, 1.0
cycles: reduce

Plotting: # Specify plots for pdf

v3: specs/ppb
O3
NO NO2
PAN CH3CHO

v3: specs/ppt
HONO
CLNO2

v3: rates/mlc
NO2-->O+NO
