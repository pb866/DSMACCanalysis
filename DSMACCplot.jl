#!/usr/local/bin/julia


"""
# Module DSMACCplot

Use DSMACC dictionary with netCDF data to plot species concentrations
and/or reaction rates for specified species/reactions and netCDF files.
"""
module DSMACCplot


##################
###  PREAMBLE  ###
##################

println("initialise...")
# Define location of external self-made modules
# (Add or modify to include your own directories)
# Local Mac:
if isdir("/Applications/bin/data/jl.mod") &&
  all(x->x!="/Applications/bin/data/jl.mod", LOAD_PATH)
  push!(LOAD_PATH,"/Applications/bin/data/jl.mod")
end
# earth0:
if isdir("~/Util/auxdata/jl.mod") &&
  all(x->x!="~/Util/auxdata/jl.mod", LOAD_PATH)
  push!(LOAD_PATH,"~/Util/auxdata/jl.mod")
end

# Assume either DSMACC main folder or DSMACC/AnalysisTools/DSMACCplot
# as current directory, other wise add/adjust folder path here:
push!(LOAD_PATH,"./jl.mod"); push!(LOAD_PATH,"AnalysisTools/DSMACCplot/jl.mod")
if splitdir(pwd())[2] == "DSMACCplot"  def_dir = "../../save/results"
else def_dir = "./save/results"
end

# Load modules/functions/python libraries
using PyPlot, DataFrames
using fhandle: test_file
using NCload
# Define the netCDF file from the first script argument
for i = 1:1-length(ARGS)  push!(ARGS,"")  end


#####################
###  MAIN SCRIPT  ###
#####################

println("load data...")
spec = []; rate = []
ifiles = split(ARGS[1])
for file in ifiles
  ncfile  = test_file(file, default_dir=def_dir)
  spc, rat = get_ncdata(ncfile)
  push!(spec,spc); push!(rate,rat)
end

println("plot data...")
plot(spec[1][:HO2],label="v3")
plot(spec[2][:HO2],label="v4")
legend()
grid()
show()
plot(rate[1][Symbol("HNO3-->NO2+OH")],label="v3")
plot(rate[2][Symbol("HNO3-->NO2+OH")],label="v4")
legend()
grid()
show()
end #module DSMACCplot
