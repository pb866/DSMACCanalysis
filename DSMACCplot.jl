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

# Assume either DSMACC main folder or DSMACC/AnalysisTools/DSMACCanalysis
# as current directory, other wise add/adjust folder path here:
push!(LOAD_PATH,"./jl.mod"); push!(LOAD_PATH,"AnalysisTools/DSMACCanalysis/jl.mod")
# push!(LOAD_PATH,".")

# Load modules/functions/python libraries
using PyPlot, DataFrames
import make_plots
using jlplot
# Define the netCDF file from the first script argument
for i = 1:2-length(ARGS)  push!(ARGS,"")  end


#####################
###  MAIN SCRIPT  ###
#####################

println("load data...")
# Read input file
commission, fidx = commission_plot(ARGS[1])
# Get nc file names and scenario names
ncfiles, label = get_scenario(commission,fidx[1])
# Read from input file, which plots to generate
icase, what, plotdata = prepare_plots(commission[fidx[2]:fidx[3]], label)
# Save DSMACC output using a python script
# to read nc files and save in Julia format
output = DSMACCoutput(ncfiles)
# Data is store in a dictionary distinguishing between "specs" and "rates"
# In each dictionary entry is is an array for each scenario (nc file)
# with dataframes for the respectives species concentrations or reaction rates


println("plot data...")
# Initilise counter for output plots
pdf = 0
# Loop over different plot types
for n = 1:length(icase)
  # Loop over different plots in each plot type
  for (i, case) in enumerate(plotdata[n])
    # Increase counter for plots and generate plots
    pdf += 1
    make_plots.lineplot(output["time"],output[what[n]],label,
                        icase[n],case,"$(lpad(pdf,4,"0")).pdf")
  end
end

# compile subplots to overall pdf
if ARGS[2] == ""  ARGS[2] = join(label,"_")  end
make_plots.compile_pdf(".",ARGS[2],ndel=pdf)

println("done.")

end #module DSMACCplot
