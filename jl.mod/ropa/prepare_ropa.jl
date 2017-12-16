__precompile__()

"""
# Module prepare_ropa

Analyse all species and reactions in a list of netCDF files and find the sinks
and sources for each species in each scenario.


## Public functions
- get_names
- split_rxn

## Private functions
- get_bratio
- merge_DUPL
"""
module prepare_ropa

  ##################
  ###  PREAMBLE  ###
  ##################

  export get_names,
         split_rxn

  using DataFrames


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    get_names(specs::Array{Any, 1}, rates::Array{Any, 1})

From the column headers in the list of dataframes for `specs` and `rates`,
get a list of species names and reactions in each secenario and return it
as list of String arrays `spc_names` and `rxns`.
"""
function get_names(specs::Array{Any, 1}, rates::Array{Any, 1})

  # Initilise overall list for species names and reactions in each scenario
  spc_names = []; rxns = []
  # Loop over scenarios
  for i = 1:length(specs)
    # Save species names and reactions for each scenario as String array
    push!(spc_names,string.(names(specs[i]))); push!(rxns,string.(names(rates[i])))
    # Exclude special cases for species
    for var in ["TIME","LAT","LON","PRESS","TEMP","H2O","M","R","RO2","EMISS","DUMMY"]
      idx = find(spc_names[i].==var);  deleteat!(spc_names[i],idx)
    end
    # Exclude special cases for reactions
    for var in ["TIME","LAT","LON","PRESS","TEMP","M","EMISS-->EMISS","R-->R"]
      idx = find(rxns[i].==var);  deleteat!(rxns[i],idx)
    end
  end

  # Return list with final species names/reactions for each scenario
  return spc_names, rxns
end #function get_names


"""
    split_rxn(rxn_array)

From a list of dataframes with the rate data in each scenario `rxn_array`,
derive an array of arrays holding matrices with branching ratios and species names
for each reaction in each scenario and return seperate arrays for educts and
products.
"""
function split_rxn(rxn_array)
  # Initilise outer array with list of educts and products for each scenario
  edarr = []; prarr = []
  # Loop over dataframes of each scenario
  for arr in rxn_array
    # Initilise inner array with list of Matrixes of branching ratios and single species
    # in each reaction
    educt = []; product = []
    # Loop over reactions in each scenario
    for rxn in arr
      # Split reactios into single educts and single products
      # (with possible branching ratios)
      edct, prdct = split.(split(rxn,"-->"),"+")
      # Get matrices with branching ratios and species names
      # for educts/products in each reaction
      edct=get_bratio(edct); prdct=get_bratio(prdct)
      edct=merge_DUPL(edct); prdct=merge_DUPL(prdct)
      # Save matrices
      push!(educt,edct); push!(product,prdct)
    end
    # Save inner array in outer array
    push!(edarr,educt); push!(prarr,product)
  end

  # Return output arrays
  return edarr, prarr
end


###########################
###  PRIVATE FUNCTIONS  ###
###########################

"""
    get_bratio(spec)

From list of branching ratio/species string `spec` derived from the netCDF reaction
list in rates, generate a matrix with branching ratios as floats and strings of
species names and return it.
"""
function get_bratio(spec)

  # Initilise matrix with species names and their branching ratios
  spc_ratio = []
  # Find branching ratios (if defined) for each ratio/species string
  ratio = match.(r"^[0-9.]+",spec)
  # Loop over ratio/species string
  for i = 1:length(spec)
    # Assign branching ratio of 1, if missing
    # and save to matrix with species name
    if ratio[i] == nothing  spc_ratio=vcat(spc_ratio,[1.0 spec[i]])
    else
      # If branching ratio is found split ratio and species and
      # save separately in matrix
      spec[i]=replace(spec[i],r"^[0-9.]+","")
      spc_ratio=vcat(spc_ratio,[float(ratio[i].match) spec[i]])
    end
  end

  # Return final matrix with branching ratio and name for each species
  return spc_ratio
end


"""
    merge_DUPL(spc_ratio)

From matrix with branching ratios and species names, find duplicate species and
merge them and their branching ratios. Return the revised matrix `spc_ratio`.
"""
function merge_DUPL(spc_ratio)

  # Loop over species
  strt = 1
  for i = length(spc_ratio[:,2]):-1:strt+1
    for j = length(spc_ratio[:,2])-1:-1:strt
      # Find duplicate species
      if spc_ratio[i,2] == spc_ratio[j,2]
        # Combine species and add up branching ratios
        spc_ratio[j,1] += spc_ratio[i,1]
        spc_ratio = spc_ratio[1:size(spc_ratio,1).!=i,:]
        strt += 1
      end
    end
  end
  # Return the adjusted mechanism data
  return spc_ratio
end #function merge_DUPL

end #module prepare_ropa
