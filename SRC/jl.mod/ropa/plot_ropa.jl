# Don't use __precompile__() with pyimport!

"""
# Module plot_ropa

Plot time-resolved sink and source analysis for a list of species in the given
scenario(s).

# Functions
## public
- plot_fluxes
- load_plotdata
- plot_data
- find_SCENidx

## private
- spc_stats
- get_Xdata
- get_Ydata
"""
module plot_ropa

##################
###  PREAMBLE  ###
##################

# Export public functions
export plot_fluxes,
       load_plotdata,
       plot_data,
       find_SCENidx

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

# Import Julia and self-made modules
using PyCall, PyPlot, DataFrames
using make_plots

# Import python functions
@pyimport numpy as np
@pyimport matplotlib.backends.backend_pdf as pdf
# np = pywrap(pyimport("numpy"))


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    plot_fluxes(species,scenarios,fname,sources,sinks,concs;
                llim::Float64=0.05, ulim::Float64=0.7)

From an array with MCM `species` names and a `scenarios` list as well as the
ropa analysis data for `sources`, `sinks`, and `concs`, plot the time-resolved
sink and source fluxes to a pdf named `fname`. An optional lower (`llim`) and
upper (`ulim`) cut-off parameter (default: 5%/70%) can be defined to combine
minor fluxes and omit major fluxes in the plots.
"""
function plot_fluxes(species,scenarios,fname,sources,sinks,concs;
                     llim::Float64=0.05, ulim::Float64=0.7)
  # Initilise pdf output file
  pdffile = pdf.PdfPages(fname)

  # Loop over species and scenarios
  for spc in species  for scen in scenarios
    # Get index for current scenario
    s = find_SCENidx(scenarios,scen)
    # Load plot data from ropa analysis
    modtime, src, snk, src_rev, snk_rev =
      load_plotdata(spc,s,sources,sinks,concs,scenarios,llim=llim,ulim=ulim)
    # Generate plots
    fig = plot_data(spc,scen,modtime,src,snk)
    if fig != nothing  pdffile[:savefig](fig)  end
    fig = plot_data(spc,scen,modtime,src_rev,snk_rev)
    if fig != nothing  pdffile[:savefig](fig)  end
  end  end
  # Close output pdf
  pdffile[:close]()
end #function plot_fluxes


"""
    load_plotdata(spc,s,sources,sinks,conc,scen;llim::Float64=0.05,ulim::Float64=0.7)

From the species MCM name `spc` and the chosen scenario with index `s` as well as the
`sources` and `sinks` flux data and the species `conc`entrations prepare plot data
for a time-resolved sink and source analysis.

A lower (`llim`) and (`upper`) cut-off (by default 5%/70%) for minor/major fluxes
can be defined. Minor fluxes are groupednin the plots, major fluxes are omitted.
Furthermore, the full list with scenario names (`scen`) is needed.
"""
function load_plotdata(spc,s,sources,sinks,conc,scen;
                       llim::Float64=0.05, ulim::Float64=0.7)
  # Generate y data for sources
  idx, fraction = spc_stats(spc,s,sources)
  if idx != nothing  Ysrc, Ysrc_rev = get_Ydata(sources,s,idx,fraction,llim,ulim)
  else
    Ysrc = ("no fluxes","no fluxes","no fluxes")
    Ysrc_rev = ("no fluxes","no fluxes","no fluxes")
  end
  # Generate y data for sinks
  idx, fraction = spc_stats(spc,s,sinks)
  if idx != nothing  Ysnk, Ysnk_rev = get_Ydata(sinks,s,idx,fraction,llim,ulim)
  else
    Ysnk = ("no fluxes","no fluxes","no fluxes")
    Ysnk_rev = ("no fluxes","no fluxes","no fluxes")
  end

  # Return model time, source and sink data (with labels as tuple)
  return Ysrc, Ysnk, Ysrc_rev, Ysnk_rev
end #function load_plotdata


"""
    plot_data(spc,scen,modtime,src,snk)

Plot sources (`src`) and sinks (`snk`) over model time as a stacked area plot
for the current species `spc` int the current scenario `scen`.
"""
function plot_data(spc,scen,modtime,src,snk,nights,pltnight,t_frmt,sfile)
  # Plot data and save plots
  if src[1]=="no fluxes" && snk[1]=="no fluxes"
    # Ignore plots with no valid data
    fig = nothing
  elseif src[1]=="no fluxes"
    # Plots without sources
    # Set sink fluxes negative
    snk[1] .*= -1.
    # Set colour scheme
    cs, dt = sel_ls(cs="sink",nc=1:length(snk[2]))
    fig = plot_flux(spc, scen, modtime, snk, cs, nights, pltnight, t_frmt, sfile)
    # Restore original sink fluxes
    snk[1] .*= -1.
  elseif snk[1]=="no fluxes"
    # Plots without sinks
    # Set colour scheme
    cs, dt = sel_ls(cs="source",nc=1:length(src[2]))
    fig = plot_flux(spc, scen, modtime, src, cs, nights, pltnight, t_frmt, sfile)
  else
    # Plots with sources and sinks
    fig = plot_prodloss(spc, scen, modtime, src, snk, nights, pltnight, t_frmt, sfile)
  end

  # Return PyObject with plot data
  return fig
end #function plot_data


"""
    find_SCENidx(scenarios, chosen_scen)

From the list of `scenarios` and the string with the chosen scenario (`chosen_scen`),
return the index of the chosen scenario used in all arrays from the ropa analysis and
for plotting.
"""
function find_SCENidx(scenarios, chosen_scen)

  # Initilise
  s = []
  # Try to get current scenario index
  try s = find(scenarios.==chosen_scen)[1]
  catch
    # Warn on failure and abort
    println("Scenario name $chosen_scen in current scenario list not found.")
    println("Script aborted."); exit()
  end
  if length(s) ≥ 2
    # Warn, if mulitple indices were found and continue with 1st one
    println("Warning! Several scenarios with the name $chosen_scen found.")
    println("First scenario chosen.")
  end

  return s
end #function find_SCENidx


###########################
###  PRIVATE FUNCTIONS  ###
###########################


"""
    spc_stats(spc,s,flux_data)

Get the index of the chosen species `spc` in the current scenario with index `s`
from the `flux_data` and return it to together with the mean fractional contribution
of each flux (either sink or source flux).
"""
function spc_stats(spc,s,flux_data)
  # Get species index
  idx=find(flux_data[s][:,1].==spc)[1]
  # Get total flux
  tot=0.
  try tot=mean(sum(flux_data[s][idx,2][:,2]))
  catch
    idx=nothing; fraction=nothing
    return idx, fraction
  end
  # Get fractional contributions
  fraction=mean.(flux_data[s][idx,2][:,2])./tot

  # Return species index and fractions
  return idx, fraction
end #function spc_stats


"""
    get_Ydata(flux_data,s,idx,fraction,llim,ulim)

From the `flux_data` (sinks or sources), the scenario index `s`, the species index
`idx`, the `fraction`s of each flux, and the cut-off parameters for minor (`llim`)
and major (`ulim`) fluxes, determine and return a tuple with the flux data for
plotting (where minor fluxes are combined) and their legend titles and a "no fluxes"
string as well as a second tuple with the revised plotting data (without major fluxes
and combined minor fluxes), the legend and a list of the omitted major fluxes.
"""
function get_Ydata(flux_data,s,idx,fraction,llim,ulim)

  # Define large fluxes (main and major)
  nmain = find(llim.≤fraction.≤ulim)
  nmax = find(fraction.>ulim)
  # Reassign data to main, if only one flux exists below ulim
  if length(fraction)==1 && ulim!=1.0
    nmain = nmax
    nmax = Int64[]
  end
  # Save main and major fluxes in arrays
  main  = flux_data[s][idx,2][nmain,:]
  major = flux_data[s][idx,2][nmax,:]
  # Assign fluxes as y data
  ydata = vcat(major,main)
  # Add accumulated minor fluxes to y data
  ydata = def_minor(ydata,flux_data[s][idx,2],fraction,llim)

  if !isempty(major)
    # If fluxes above upper limit exist, exclude them and recalculate fractions
    rev_fluxes = vcat(main,flux_data[s][idx,2][find(fraction.<llim),:])
    rev_fraction=mean.(rev_fluxes[:,2])./sum(mean.(rev_fluxes[:,2]))
    # Define main fluxes as y data
    rev_nmain = find(rev_fraction.≥llim)
    yrev = rev_fluxes[rev_nmain,:]
    # Add accumulated minor fluxes to y data
    yrev = def_minor(yrev,rev_fluxes,rev_fraction,llim)
    yrev_data = np.row_stack(tuple(yrev[:,2]))
  else
    # Set yref to nothing, if no fluxes exist above upper threshold
    yrev_data = "no fluxes"; yrev = ["no fluxes"]; major = ["no fluxes"]
  end

  # Return adjusted fluxes and their legends as a tuple
  return (np.row_stack(tuple(ydata[:,2])),ydata[:,1],"no fluxes"), (yrev_data, yrev[:,1], major[:,1])
end #function get_Ydata


"""
    def_minor(ydata,fluxes,fraction,llim)

Append `ydata` with main and major fluxes by combined minor fluxes using the data
in `fluxes`, the mean `fraction` of each flux and the cut-off `llim` below which
fluxes are accumulated.
"""
function def_minor(ydata,fluxes,fraction,llim)
  # Find fluxes below lower threshold
  nmin = find(fraction.<llim)
  minor = fluxes[nmin,:]
  # Combine minor fluxes and add to y data (if present)
  if !isempty(minor)
    minor = sum(minor[:,2])
    ydata = vcat(ydata, ["minor fluxes" [minor]])
  end

  # Return appended y data
  return ydata
end #function def_minor

end #module plot_ropa
