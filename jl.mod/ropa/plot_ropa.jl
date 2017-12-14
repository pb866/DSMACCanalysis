# Don't use __precompile__() with pyimport!

"""
"""
module plot_ropa

# Export public functions
export load_plotdata, plot_fluxes

# Import Julia and self-made modules
using PyCall, PyPlot, DataFrames
using make_plots

# Import python functions
# @pyimport numpy as np
@pyimport matplotlib.backends.backend_pdf as pdf
np = pywrap(pyimport("numpy"))


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    plot_fluxes(species,scenarios,fname,sources,sinks,concs,;cut_off::Float64=0.05)

From an array with MCM `species` names and a `scenarios` as well as the ropa
analysis data for `sources`, `sinks`, and `concs`, plot the time-resolved sink
and source fluxes to a pdf named `fname`. An optional `cut_off` (default: 5%)
can be defined to combine minor fluxes in the plots.
"""
function plot_fluxes(species,scenarios,fname,sources,sinks,concs;
                     cut_off::Float64=0.05)
  # Initilise pdf output file
  pdffile = pdf.PdfPages(fname)

  # Loop over species and scenarios
  for spc in species  for scen in scenarios
    # Load plot data from ropa analysis
    modtime, src, snk = load_plotdata(spc,scen,sources,sinks,concs,scenarios;cut_off=0.05)
    # Plot data and save plots
    fig = plot_flux(modtime, src, snk)
    pdffile[:savefig](fig)
  end  end
  # Close output pdf
  pdffile[:close]()
end


"""
    load_plotdata(spc,chosen_scen,sources,sinks,conc,scen;cut_off::Float64=0.05)

From the species MCM name `spc` and the scenario name `chosen_scen` as well as the
`sources` and `sinks` flux data and the species `conc`entrations prepare plot data
for a time-resolved sink and source analysis.

A cut-off (by default 5%) for minor fluxes can be defined. Minor fluxes are grouped
in the plots. Furthermore, the full list with scenario names is needed.
"""
function load_plotdata(spc,chosen_scen,sources,sinks,conc,scen;
                       cut_off::Float64=0.05)
  # Get index for current scenario
  s = find_SCENidx(scen,chosen_scen)
  # Generate x and y data for sinks and sources
  modtime = get_Xdata(conc,s)
  idx, fraction = spc_stats(spc,s,sources)
  Ysrc = get_Ydata(sources,s,idx,fraction,cut_off)
  idx, fraction = spc_stats(spc,s,sinks)
  Ysnk = get_Ydata(sinks,s,idx,fraction,cut_off)
  # Negatate sink data
  Ysnk[1] .*= -1.

  # Return model time, source and sink data (with labels as tuple)
  return modtime, Ysrc, Ysnk
end


###########################
###  PRIVATE FUNCTIONS  ###
###########################

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
end


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
  tot=mean(sum(flux_data[s][idx,2][:,2]))
  # Get fractional contributions
  fraction=mean.(flux_data[s][idx,2][:,2])./tot

  # Return species index and fractions
  return idx, fraction
end #function spc_stats


"""
    get_Xdata(conc,s)

From the species concentration array `conc` and the scenario index `s`, calculate
the model time (x values) for the current plots.
"""
function get_Xdata(conc,s)
  # Initilise
  modtime = Float64[]; t = 0.
  # Get time step from last 2 times
  dt = conc[s][:TIME][end] - conc[s][:TIME][end-1]
  # Use time step to fill the array for the length of the model run with times
  # starting at 0
  for i = 1:length(conc[s][:TIME])
    push!(modtime,t); t += dt/3600
  end

  # Return time array
  return modtime
end


"""
    get_Ydata(flux_data,s,idx,fraction,cut_off)

From the `flux_data` (sinks or sources), the scenario index `s`, the species index
`idx`, the `fraction`s of each flux, and the `cut_off` parameter for minor fluxes,
determine and return a tuple with the flux data for plotting (where minor fluxes
are combined) and their legend titles.
"""
function get_Ydata(flux_data,s,idx,fraction,cut_off)

  # Find and combine minor fluxes
  excl=find(x->x<cut_off,fraction)
  rest = sum(flux_data[s][idx,2][excl,2])
  # Find major fluxes
  major = flux_data[s][idx,2][find(x->x≥cut_off,fraction),:]

  # Store both fluxes (and their legends)
  ydata = vcat(major,["minor fluxes" [rest]])

  # Return adjusted fluxes and their legends as a tuple
  return (np.row_stack(tuple(ydata[:,2])),ydata[:,1])
end

end #module plot_ropa
