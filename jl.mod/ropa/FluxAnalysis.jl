__precompile__()

"""
# Module FluxAnalysis

Determine chemical source and sink fluxes for each species in the DSMACC model
output.

# Functions
## Public functions
- flux_rates
- net_flux
- prodloss

## Private functions
- get_fluxes
"""
module FluxAnalysis

export flux_rates,
       net_flux,
       prodloss


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    flux_rates(ed_scen,specs,rates,rxn)

From `rates` with rate constants for each reaction in each scenario, the labels
for each reaction in each scenario `rxn`, the concentrations for each species in
each scenario `specs`, and the educts with their respective branching ratios in
each scenario `ed_scen` determine the chemical fluxes for each reaction in each
scenario `flux_data` and return it.
"""
function flux_rates(ed_scen,specs,rates,rxn)

  # Initilise complete flux_data array
  flux_data = []
  # Loop over scenarios
  for i = 1:length(ed_scen)
    # Initilise flux data for current scenario
    flux_rate = Matrix{Any}(0,2)
    # Loop over reactions in each scenario
    for j = 1:length(ed_scen[i])
      # Initilise chemical flux as v = k[first educt]
      v = rates[i][Symbol(rxn[i][j])].*specs[i][Symbol(ed_scen[i][j][1,2])].^ed_scen[i][j][1,1]
      # Loop over further educts in each reaction
      for k = 2:length(ed_scen[i][j][:,2])
        # Add concentrations to the product to complete calculation of chemical flux
        v .*= specs[i][Symbol(ed_scen[i][j][k,2])].^ed_scen[i][j][k,1]
      end
      # Save flux of current reaction in the flux array for current scenario
      flux_rate = vcat(flux_rate, [rxn[i][j] [v]])
    end
    # Save current scenario in overall flux data array
    push!(flux_data,flux_rate)
  end

  # Return
  return flux_data
end


"""
    net_flux(educt,product,flux_data)

From the `educt` and `product` array with species names and their branching ratios
in each scenario, determine equilibria in the `flux_data` and adjust it to use
net fluxes (negative values will be set to zero).
"""
function net_flux(educt,product,flux_data)

  # Loop over scenarios
  for s = 1:length(educt)
    # Loop over LHS and RHS of reactions and find equilibria
    # (where educts = products)
    for i=1:length(educt[s])-1
      for j=i+1:length(product[s])
        # Find equilibria
        if sort(educt[s][i][:,2])==sort(product[s][j][:,2]) && sort(product[s][i][:,2])==sort(educt[s][j][:,2])
          # Calculate net_flux by substracting flux data of LHS from RHS
          net = flux_data[s][i,2].-flux_data[s][j,2]
          # Initilise revised flux data
          flux_data[s][i,2] = zeros(length(net)); flux_data[s][j,2] = zeros(length(net))
          # Find positive values for forward and backward reaction
          f=find(net.>0); b=find(net.<0)
          # Save positive values to revised net fluxes leaving only the net production
          # or loss for the respective reaction
          flux_data[s][i,2][f] = net[f]; flux_data[s][j,2][b] = net[b].⋅-1
          # Replace arrow signs in net reactions to indicate the data manipulation
          flux_data[s][i,1]=replace(flux_data[s][i,1],"-->","<-->>")
          flux_data[s][j,1]=replace(flux_data[s][j,1],"-->","<-->>")
        end
      end
    end
  end

  # Return revised flux data
  return flux_data
end


"""
    prodloss(spc_arr,edct,prdct,flux_data)

From the species names in each scenario `spc_arr`, the species names and branching
ratios in each reaction and scenario for educts (`edct`) and products (`prdct`),
respectively, and the `flux_data`, define arrays for `sources` and `sinks` for every
species in every scenario.
"""
function prodloss(spc_arr,edct,prdct,flux_data)

  sources = get_fluxes(prdct,edct,flux_data,spc_arr)
  sinks   = get_fluxes(edct,prdct,flux_data,spc_arr)

  return sources, sinks
end


###########################
###  PRIVATE FUNCTIONS  ###
###########################

"""
    get_fluxes(reactant,other_reactant,rates,spc)

Depending on whether the `reactant` is educt or product and with the help of the
`other_reactant` (product/educt) as well as the flux data `rates` and the species
names `spc` in each scenario, return an array `fluxes_scen` with either sinks or
sources for every species in every scenario.
"""
function get_fluxes(reactant,other_reactant,rates,spc)

  # Initilise array for all turnovers in each scenario
  fluxes_scen = []
  # Loop over scenarios and initialise flux arrays for each scenario
  for s = 1:length(rates)
    fluxes = Matrix{Any}(length(spc[s]),2)
    for i = 1:length(spc[s])
      fluxes[i,:] = [spc[s][i] [Matrix{Any}(0,2)]]
    end
    # Loop over reactions in each scenario
    for i = 1:length(rates[s][:,1])
      # Loop over reactants in each reactions (either LHS or RHS)
      for j = 1:length(reactant[s][i][:,2])
        # Find reactant in species array (exclude special cases such as DUMMY)
        try idx = find(spc[s].==reactant[s][i][j,2])[1]
          # Check for the same reactant on the other side of the equation
          try ind = find(reactant[s][i][j,2].==other_reactant[s][i][:,2])[1]
            # Calculate net stoichiometric index, if reactant is on both sides
            net_br = reactant[s][i][j,1] - other_reactant[s][i][ind,1]
            # Save fluxes for positive stoichiometric indices
            if net_br ≤ 0  continue
            else fluxes[idx,2] =
                    vcat(fluxes[idx,2],[rates[s][i,1] [rates[s][i,2].⋅net_br]])
            end
          catch
            # Save fluxes modified by stoichiometric indices for general cases
            fluxes[idx,2] = vcat(fluxes[idx,2],
              [rates[s][i,1] [rates[s][i,2].⋅reactant[s][i][j,1]]])
          end
        end
      end
    end
    # Save fluxes for current scenario to overall array
    push!(fluxes_scen,fluxes)
  end

  # Return an array with sources or sinks depending on reactant for in each scenario
  return fluxes_scen
end

end #module FluxAnalysis
