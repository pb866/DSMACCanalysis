#!/usr/local/bin/julia


"""
# Module ROPAplot

Plot sink and source fluxes from ROPA analysis as bar plots for mean values and
as time-resolved sink and source analysis.
"""
module ROPAplot


##################
###  PREAMBLE  ###
##################

println("initialise...")
# Define location of external self-made modules
# (Add or modify to include your own directories)
# Local Mac:
if isdir("/Applications/bin/data/jl.mod") &&
  all(LOAD_PATH.!="/Applications/bin/data/jl.mod")
  push!(LOAD_PATH,"/Applications/bin/data/jl.mod")
end
# earth0:
if isdir("~/Util/auxdata/jl.mod") && all(LOAD_PATH.!="~/Util/auxdata/jl.mod")
  push!(LOAD_PATH,"~/Util/auxdata/jl.mod")
end

# Add path of internal self-made modules
if isdir(joinpath(Base.source_dir(),"jl.mod")) &&
  all(LOAD_PATH.!=joinpath(Base.source_dir(),"jl.mod"))
  push!(LOAD_PATH,joinpath(Base.source_dir(),"jl.mod"))
end
if isdir(joinpath(Base.source_dir(),"jl.mod/ropa")) &&
  all(LOAD_PATH.!=joinpath(Base.source_dir(),"jl.mod/ropa"))
  push!(LOAD_PATH,joinpath(Base.source_dir(),"jl.mod/ropa"))
end

# Import Julia and self-made modules
using DataFrames#, PyCall, PyPlot
using ropa_tool#, make_plots
using plot_ropa

# Import python functions
# @pyimport matplotlib.backends.backend_pdf as pdf


# Define the netCDF files from the first script argument
for i = 1:2-length(ARGS)  push!(ARGS,"")  end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

#####################
###  MAIN SCRIPT  ###
#####################

# Perform ROPA analysis on scenarios given in script argument 1
sources, sinks, concs = ropa(ARGS[1])
# Retrieve scenario names as list of file base names
scenarios = String[]
for s in basename.(split(ARGS[1]))  push!(scenarios,splitext(s)[1])  end
# Define a list of species to be plotted
species = ["HO2", "OH", "NO", "NO2", "O3"]

# Define output file name
if ARGS[2] == ""  ARGS[2] = "fluxes_"*join(scenarios,"_")  end
fname = ARGS[2]*".pdf"

# Generate file with plots form species and scenario list
plot_fluxes(species,scenarios,fname,sources,sinks,concs)

end #module ROPAplot
