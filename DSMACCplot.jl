#!/usr/local/bin/julia


"""
# Module DSMACCplot

Use DSMACC dictionary with netCDF data to plot species concentrations
and/or reaction rates for specified species/reactions and netCDF files.

If input file asks for flux plots, perform ROPA analysis and print time-resolved
sink and source analysis plots.

Group species according to the following properties, which can be outputted in
stacked area or line plots:
- Molecule's size
- Chromophore class
- O:C ratio
"""
module DSMACCplot

# Preamble with module loading
include("SRC/preamble.jl")
# Routines to load and analyse data
include("SRC/load.jl")
# Loop over species in input file and produce plots
include("SRC/plot.jl")

println("done.")

end #module DSMACCplot
