#!/usr/local/bin/julia


"""
# Module ropa

Analyse sink and source fluxes from DSMACC output.
"""
module ropa


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

# Load modules/functions/python libraries
using PyPlot, DataFrames
using NCload
using prepare_ropa, FluxAnalysis

# Define the netCDF files from the first script argument
for i = 1:1-length(ARGS)  push!(ARGS,"")  end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#


#####################
###  MAIN SCRIPT  ###
#####################

println("load data...")
# Load DSMACC output
specs = []; rates = []
for file in split(ARGS[1])
  spec, rate = NCload.get_ncdata(file)
  push!(specs,spec); push!(rates,rate)
end

println("analyse reactions...")
# Get relevant species and reactions for each scenario as array of arrays
spc, rxn = get_names(specs,rates)
# Find species/branching ratios of LHS/RHS of reaction
educt, product = split_rxn(rxn)

println("determine fluxes...")
# Calculate turnovers using net fluxes for equilibria
flux_data = flux_rates(educt,specs,rates,rxn)
flux_data = net_flux(educt,product,flux_data)
# Find sinks and sources for each species
source, sink = prodloss(spc,educt,product,flux_data)

end #module ropa
