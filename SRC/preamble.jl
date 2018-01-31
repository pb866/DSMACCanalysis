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
if isdir(joinpath(homedir(),"Util/auxdata/jl.mod")) &&
  all(LOAD_PATH.!=joinpath(homedir(),"Util/auxdata/jl.mod"))
  push!(LOAD_PATH,joinpath(homedir(),"Util/auxdata/jl.mod"))
end

# Add path of internal self-made modules
mdir = joinpath(Base.source_dir(),"jl.mod")
rdir = joinpath(Base.source_dir(),"jl.mod/ropa")
if all(LOAD_PATH.!=mdir)  push!(LOAD_PATH,mdir)  end
if all(LOAD_PATH.!=rdir)  push!(LOAD_PATH,rdir)  end

# Load modules/functions/python libraries
using PyPlot, PyCall, DataFrames
using Juno: input
import make_plots
using groupSPC
using jlplot
using ropa_tool, plot_ropa
using NCload

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

# Load python function for multiple plots in one pdf
@pyimport matplotlib.backends.backend_pdf as pdf

# Define the netCDF file(s) from the first script argument
# For several files, use list of whitespace separated files wrapped in double quotes
for i = 1:2-length(ARGS)  push!(ARGS,"")  end

# Assume plot input file in parent folder to repo, if no directory is specified
if dirname(ARGS[1]) == ""
  ARGS[1] = normpath(joinpath(Base.source_dir(),"../..",ARGS[1]))
end
println(normpath(joinpath(Base.source_dir(),"../..",ARGS[1])))
