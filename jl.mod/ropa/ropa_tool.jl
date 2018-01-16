__precompile__()

"""
# Module ropa_tool

Analyse sink and source fluxes from DSMACC output.
"""
module ropa_tool


##################
###  PREAMBLE  ###
##################

export ropa

# Find directory of module source code
cdir=Base.source_dir()
# Find parent directory
pdir=dirname(cdir)

# Define location of external self-made modules
# (Add or modify to include your own directories)if all(LOAD_PATH.!=cdir)  push!(LOAD_PATH,cdir)  end
if all(LOAD_PATH.!=cdir)  push!(LOAD_PATH,cdir)  end
if all(LOAD_PATH.!=pdir)  push!(LOAD_PATH,pdir)  end

# Load modules/functions/python libraries
using DataFrames
using NCload
using prepare_ropa, FluxAnalysis

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#


###################
###  FUNCTIONS  ###
###################

"""
    function ropa(scenarios=[]; specs=[], rates=[], cycles="reduce")

Determine species concentrations `specs` and rate constants `rates`, if not provided
as function parameters and return the sources and sinks as array for each species
in each scenario as well as the species concentrations.

Reduce inorganic NOx and Ox cycles to main sources and sinks, if not specified
otherwise by parameter `cycles`.
"""
function ropa(scenarios=[]; specs=[], rates=[], cycles="reduce")
  println("load ROPA data...")
  # Load DSMACC output
  if isempty(specs) || isempty(rates)
    for file in scenarios
      spec, rate = get_ncdata(file)
      push!(specs,spec); push!(rates,rate)
    end
  end

  println("analyse reactions...")
  # Get relevant species and reactions for each scenario as array of arrays
  spc, rxn = get_names(specs,rates)
  # Find species/branching ratios of LHS/RHS of reaction
  educt, product = split_rxn(rxn)

  println("determine fluxes...")
  # Calculate turnovers using net fluxes for equilibria
  flux_data = flux_rates(educt,specs,rates,rxn)
  # Calculate net fluxes from cycles and equilibria
  if cycles=="reduce"
    educt,product,flux_data = net_cycles(educt,product,flux_data)
  end
  flux_data = net_flux(educt,product,flux_data)
  # Delete fluxes with negligible turnovers
  educt,product,flux_data = del_zerofluxes(educt,product,flux_data)
  # Find sinks and sources for each species
  source, sink = prodloss(spc,educt,product,flux_data)

  # Return sink and source arrays
  return source, sink, specs
end #function ropa

end #module ropa_tool
