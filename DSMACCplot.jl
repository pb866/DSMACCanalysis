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
using PyPlot, PyCall, DataFrames
import make_plots
using groupSPC
using jlplot
# Define the netCDF file from the first script argument
for i = 1:2-length(ARGS)  push!(ARGS,"")  end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

#####################
###  MAIN SCRIPT  ###
#####################

println("load data...")
# Read input file
commission, fidx = commission_plot(ARGS[1])
# Get nc file names and scenario names
ncfiles, label = get_scenario(commission,fidx[1])
# Read from input file, which plots to generate
icase, what, unit, plotdata = prepare_plots(commission[fidx[2]:fidx[3]], label)
# Save DSMACC output using a python script
# to read nc files and save in Julia format
output = DSMACCoutput(ncfiles)
# Data is store in a dictionary distinguishing between "specs" and "rates"
# In each dictionary entry is is an array for each scenario (nc file)
# with dataframes for the respectives species concentrations or reaction rates

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

println("analyse data...")
# Load database with different species names
spcDB = readDB("IO/MCMv331species.db")
# Translate species names from MCM to GECKO-A
gspc = translateNMVOC(output["specs"],spcDB)
# Group species by properties
CC, OC, CN, chrom_class, OCratio_class, size_class = group_specs(gspc,spcDB)
# Add new classes to dataframes
output["specs"] = add_conc(output["specs"],chrom_class,OCratio_class,size_class,CC,OC,CN)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

println("plot data...")
# Initilise counter for output plots
if ARGS[2] == ""  ARGS[2] = join(label,"_")*".pdf"  end
@pyimport matplotlib.backends.backend_pdf as pdf
figures = []
# Loop over different plot types
for n = 1:length(icase)
  # Loop over different plots in each plot type
  for (i, case) in enumerate(plotdata[n])
    # Increase counter for plots and generate plots
    fig = make_plots.lineplot(output["time"],output[what[n]],label,what[n],unit[n],
                              icase[n],case)
    push!(figures,fig)
  end
end
pdffile = pdf.PdfPages(ARGS[2]) # create pdf file
[pdffile[:savefig](f) for f in figures] # add figures to file
pdffile[:close]() # close pdf file


println("done.")

end #module DSMACCplot
