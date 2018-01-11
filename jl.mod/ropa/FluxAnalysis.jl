__precompile__()

"""
# Module FluxAnalysis

Determine chemical source and sink fluxes for each species in the DSMACC model
output. If parameter `cycles` is set to `"reduce"`, calculate net cycles for the
inorganic NOx and Ox cycles.

# Functions
## Public functions
- flux_rates
- net_cycles
- net_flux
- del_zerofluxes
- prodloss

## Private functions
- get_fluxes
- init_maincycles
- calc_maincycles
- del_negval
- del_original
"""
module FluxAnalysis

export flux_rates,
       net_cycles,
       net_flux,
       del_zerofluxes,
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
end #function flux_rates


"""
    net_cycles(educt,product,flux_data)

From the `flux_data` and the respective `educt`/`product` lists, determine net
fluxes for the major inorganic NOx and Ox cycles for O3, O1D, O3P, NO, and NO2
and adjust the `flux_data` array by adding new net fluxes and deleting the
individual fluxes. `educt` and `product` arrays are adjusted accordingly.
"""
function net_cycles(educt,product,flux_data)

  # Add array entries in flux_data for net fluxes
  educt,product,flux_data,fluxlen=init_maincycles(educt,product,flux_data)
  # Calculate net fluxes and delete individual fluxes
  # and the corresponding educts/products
  educt,product,flux_data=calc_maincycles(educt,product,flux_data,fluxlen)

  # Return adjusted arrays
  return educt,product,flux_data
end #function net_cycles


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
end #function net_flux


"""
    del_zerofluxes(flux_data)

Delete fluxes in `flux_data` with zero turnover.
"""
function del_zerofluxes(educt,product,flux_data)
  # Loop over scenarios
  for s = 1:length(flux_data)
    # Find and keep fluxes with non-negligible turnovers
    keep=find(sum.(flux_data[s][:,2]).≥1.e-10)
    flux_data[s] = flux_data[s][keep,:]
    # Additionally, delete educts and products, if a flux is 0
    educt[s]     = educt[s][keep]
    product[s]   = product[s][keep]
  end

  # Return revised flux_data
  return educt,product,flux_data
end #function del_zerofluxes


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
end #function prodloss


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
      # Loop over reactants in each reaction (either LHS or RHS)
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
end #function get_fluxes


"""
    init_maincycles(educt,product,flux_data)

Add empty entries for net cycles to the `flux_data` array and add the respective
source and sink species to the `educt` and `product` lists.
"""
function init_maincycles(educt,product,flux_data)
  # Initilise original vector lengths for each scenario
  fluxlen=[]
  for s=1:length(flux_data)
    # Save vector length for current scenario
    push!(fluxlen,length(flux_data[s][:,1]))
    # Initilise net fluxes of main inorganic cycles
    flux_data[s]=vcat(flux_data[s],["O3 net main prod" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["O3 net main loss" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["O net main prod" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["O net main loss" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["O1D net main prod" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["O1D net main loss" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["NO net main prod" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["NO net main loss" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["NO2 net main prod" [zeros(flux_data[s][1,2])]])
    flux_data[s]=vcat(flux_data[s],["NO2 net main loss" [zeros(flux_data[s][1,2])]])

    # Initilise educt/product arrays for main net cycles
    push!(educt[s],[1.0 ""]); push!(product[s],[1.0 "O3"])
    push!(educt[s],[1.0 "O3"]); push!(product[s],[1.0 ""])
    push!(educt[s],[1.0 ""]); push!(product[s],[1.0 "O"])
    push!(educt[s],[1.0 "O"]); push!(product[s],[1.0 ""])
    push!(educt[s],[1.0 ""]); push!(product[s],[1.0 "O1D"])
    push!(educt[s],[1.0 "O1D"]); push!(product[s],[1.0 ""])
    push!(educt[s],[1.0 ""]); push!(product[s],[1.0 "NO"])
    push!(educt[s],[1.0 "NO"]); push!(product[s],[1.0 ""])
    push!(educt[s],[1.0 ""]); push!(product[s],[1.0 "NO2"])
    push!(educt[s],[1.0 "NO2"]); push!(product[s],[1.0 ""])
  end

  # Return ammended arrays for reactants and fluxes and original array lengths
  return educt,product,flux_data,fluxlen
end #function init_maincycles


"""
    calc_maincycles(educt,product,flux_data,fluxlen)

Calculate net fluxes for inorganic NOx and Ox species in the flux_data array
up to the index `fluxlen` (holding the individual fluxes without the net fluxes)
using the `educt` and `product` lists to identify the individual fluxes and
equilibria.
"""
function calc_maincycles(educt,product,flux_data,fluxlen)
  # Loop over scenarios
  for s = 1:length(flux_data)
    # Loop over fluxes in reverse order (to be able to delete fluxes)
    # starting at the individual fluxes without the initialised net flux part
    for i = fluxlen[s]:-1:1
      # Calculate net fluxes for major ozone, O(1D), O(3P), NO and NO3 sources and sinks
      if all(educt[s][i].==[1.0 "O3"]) && all(product[s][i].==[1.0 "O1D"])
        # Add fluxes of reaction O3 -> O(1D) to the respective net cycles
        flux_data[s][end-9,2]=flux_data[s][end-9,2]-flux_data[s][i,2]
        flux_data[s][end-8,2]=flux_data[s][end-8,2]+flux_data[s][i,2]
        flux_data[s][end-5,2]=flux_data[s][end-5,2]+flux_data[s][i,2]
        flux_data[s][end-4,2]=flux_data[s][end-4,2]-flux_data[s][i,2]
        # Delete individual fluxes and the respective educt/product entries
        educt,product,flux_data=del_original(s,i,educt,product,flux_data)
      elseif all(educt[s][i].==[1.0 "O1D"]) && all(product[s][i].==[1.0 "O"])
        # Add fluxes of reaction O(1D) -> O(3P) to the respective net cycles
        flux_data[s][end-7,2]=flux_data[s][end-7,2]+flux_data[s][i,2]
        flux_data[s][end-6,2]=flux_data[s][end-6,2]-flux_data[s][i,2]
        flux_data[s][end-5,2]=flux_data[s][end-5,2]-flux_data[s][i,2]
        flux_data[s][end-4,2]=flux_data[s][end-4,2]+flux_data[s][i,2]
        # Delete individual fluxes and the respective educt/product entries
        educt,product,flux_data=del_original(s,i,educt,product,flux_data)
      elseif all(educt[s][i].==[1.0 "O"]) && all(product[s][i].==[1.0 "O3"])
        # Add fluxes of reaction O(3P) -> O3 to the respective net cycles
        flux_data[s][end-9,2]=flux_data[s][end-9,2]+flux_data[s][i,2]
        flux_data[s][end-8,2]=flux_data[s][end-8,2]-flux_data[s][i,2]
        flux_data[s][end-7,2]=flux_data[s][end-7,2]-flux_data[s][i,2]
        flux_data[s][end-6,2]=flux_data[s][end-6,2]+flux_data[s][i,2]
        # Delete individual fluxes and the respective educt/product entries
        educt,product,flux_data=del_original(s,i,educt,product,flux_data)
      elseif all(educt[s][i].==[1.0 "O3"]) && all(product[s][i].==[1.0 "O"])
        # Add fluxes of reaction O3 -> O(3P) to the respective net cycles
        flux_data[s][end-9,2]=flux_data[s][end-9,2]-flux_data[s][i,2]
        flux_data[s][end-8,2]=flux_data[s][end-8,2]+flux_data[s][i,2]
        flux_data[s][end-7,2]=flux_data[s][end-7,2]+flux_data[s][i,2]
        flux_data[s][end-6,2]=flux_data[s][end-6,2]-flux_data[s][i,2]
        # Delete individual fluxes and the respective educt/product entries
        educt,product,flux_data=del_original(s,i,educt,product,flux_data)
      elseif all(educt[s][i].==[1.0 "NO2"]) &&
        all(any(hcat(all(product[s][i].==[1.0 "O"],2),all(product[s][i].==[1.0 "NO"],2)),1))
        # Add fluxes of reaction NO2 -> NO + O(3P) to the respective net cycles
        flux_data[s][end-7,2]=flux_data[s][end-7,2]+flux_data[s][i,2]
        flux_data[s][end-6,2]=flux_data[s][end-6,2]-flux_data[s][i,2]
        flux_data[s][end-3,2]=flux_data[s][end-3,2]+flux_data[s][i,2]
        flux_data[s][end-2,2]=flux_data[s][end-2,2]-flux_data[s][i,2]
        flux_data[s][end-1,2]=flux_data[s][end-1,2]-flux_data[s][i,2]
        flux_data[s][end,2]=flux_data[s][end,2]+flux_data[s][i,2]
        # Delete individual fluxes and the respective educt/product entries
        educt,product,flux_data=del_original(s,i,educt,product,flux_data)
      elseif all(any(hcat(all(educt[s][i].==[1.0 "O3"],2),all(educt[s][i].==[1.0 "NO"],2)),1)) &&
        all(product[s][i].==[1.0 "NO2"])
        # Add fluxes of reaction NO + O3 -> NO2 to the respective net cycles
        flux_data[s][end-9,2]=flux_data[s][end-9,2]-flux_data[s][i,2]
        flux_data[s][end-8,2]=flux_data[s][end-8,2]+flux_data[s][i,2]
        flux_data[s][end-3,2]=flux_data[s][end-3,2]-flux_data[s][i,2]
        flux_data[s][end-2,2]=flux_data[s][end-2,2]+flux_data[s][i,2]
        flux_data[s][end-1,2]=flux_data[s][end-1,2]+flux_data[s][i,2]
        flux_data[s][end,2]=flux_data[s][end,2]-flux_data[s][i,2]
        # Delete individual fluxes and the respective educt/product entries
        educt,product,flux_data=del_original(s,i,educt,product,flux_data)
      end
    end #i: fluxes
    # Delete negative values in net cylcles
    flux_data[s] = del_negval(flux_data[s])
  end #s: scenarios
  # Return ammended arrays for reactants and fluxes
  return educt,product,flux_data
end #function calc_maincycles


"""
    del_negval(flux_data)

Delete negative values in the `flux_data` of net reactions.
"""
function del_negval(flux_data)
  # Find length of array in current scenario and loop over last 10 entries
  ln = length(flux_data[:,1])
  for i = ln-9:ln
    # Find negative values in net cycles
    neg = find(flux_data[i,2].<0)
    # Delete negative values in net cycles
    flux_data[i,2][neg] = 0.0
  end

  # Return adjusted flux data array
  return flux_data
end #function del_negval


"""
    del_original(s,n,educt,product,flux_data)

Delete individual fluxes with index `n` in `flux_data`, when net fluxes where
calculated, and the corresponding entries for `educt`s and `product`s for the
scenario with index `s` and return the adjusted arrays for `flux_data`, `educt`,
and `product`.
"""
function del_original(s,n,educt,product,flux_data)

  # Delete current row in flux_data matrix
  keep = find(collect(1:length(flux_data[s][:,1])).≠n)
  flux_data[s] = flux_data[s][keep,:]
  # Delete current row in reactants matrices
  deleteat!(educt[s],n); deleteat!(product[s],n)

  # Return adjusted matrices
  return educt,product,flux_data
end #function del_original

end #module FluxAnalysis
