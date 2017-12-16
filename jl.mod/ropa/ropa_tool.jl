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

# Load modules/functions/python libraries
using DataFrames
using NCload
using prepare_ropa, FluxAnalysis

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#


###################
###  FUNCTIONS  ###
###################

function ropa(scenarios)
  println("load data...")
  # Load DSMACC output
  specs = []; rates = []
  for file in split(scenarios)
    spec, rate = get_ncdata(file)
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

  # Return sink and source arrays
  return source, sink, specs
end #function ropa

end #module ropa_tool
