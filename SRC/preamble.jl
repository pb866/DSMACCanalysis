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
# Routines linked to python
try using PyPlot
catch
  Pkg.add("PyPlot")
  using PyPlot
end
using PyCall
# Data handling and data input
try using DataFrames
catch
  Pkg.add("DataFrames")
  using DataFrames
end
try using Juno: input
catch
  Pkg.add("Juno")
  using Juno: input
end
using NCload
# Plotting
try import make_plots
catch
  download("https://raw.githubusercontent.com/pb866/auxdata/master/jl.mod/make_plots.jl",
           joinpath(Base.source_dir(),"jl.mod","make_plots.jl"))
  import make_plots
end
using jlplot
# Data analysis
using groupSPC
using ropa_tool, plot_ropa

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
