Scenarios: # Specify netCDF files and scenario labels
DATA/TOLm3.nc DATA/TOLm4.nc
"MCMv3.3.1" "MCM/GECKO-A"

# General plotting parameters
Settings:
MCM:     v3.3.1
time:    TIME
night:   black/0.2
cut-off: 0.05, 1.0
cycles:  reduce
Fig:     on


Plotting: # Specify plots
# Concentration of primary compound
MCMv3.3.1, MCM/GECKO-A: specs/ppb
TOLUENE
GLYOX

# Important unsaturated dicarbonyl intermediates
MCMv3.3.1, MCM/GECKO-A: specs/ppt
MALDIAL, C5DICARB

#Oxidant concentrations
MCMv3.3.1, MCM/GECKO-A: specs/ppb
O3
NO NO2

MCMv3.3.1, MCM/GECKO-A: specs
OH
HO2

# Stacked area plots of lumped concentrations
MCMv3.3.1 MCM/GECKO-A: stack/ppb
Ald Ket DiCar Kete Nitro DiNitro Nit DiNit PNit DiPNit PAN ROOH DiROOH PAA SCI Poly NoChr Rad New
oc1 oc2 oc3 oc4
sc1 sc2 sc3 sc4 sc5 sc6 sc7 sc8 sc9 sc0

# Source and sink fluxes from ROPA analysis
MCMv3.3.1, MCM/GECKO-A: fluxes
TOLUENE
OH, HO2

MCMv3.3.1, MCM/GECKO-A: fluxes/ppb
O3
NO, NO2

# Reaction rates
MCM/GECKO-A: rates
C4MDIAL-->0.5newCYC1+0.5PXYFUONE
PR1O2HNO3-->0.5newRAD+0.5PRONO3AO+0.5NO2+0.5OH
TOLUENE+OH-->C6H5CH2O2  TOLUENE+OH-->CRESOL+HO2

Comments:
MALDIAL, C5DICARB
GLYOX
NO, NO2
TLOBIPEROH

# Stacked area plots of lumped concentrations
MCMv3.3.1 MCM/GECKO-A: stack
Ald Ket DiCar Kete Nitro DiNitro Nit DiNit PNit DiPNit PAN ROOH DiROOH PAA SCI Poly NoChr Rad New
oc1 oc2 oc3 oc4
sc1 sc2 sc3 sc4 sc5 sc6 sc7 sc8 sc9 sc0


Comments: # Any additional comments here


Above script shows example plots for a numerical experiment oxidising
10ppb of TOLUENE using 40ppb of ozone and 10ppb of NO.

Further script developments:
• ropa (bar plot of mean sources/sinks from ropa analysis)
• Unit conversions for flux plots to ppm/h, ppb/h, and ppt/h
• Grey shades for night-time hours (needs calculation of sza)
