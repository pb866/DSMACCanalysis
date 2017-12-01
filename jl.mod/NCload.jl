#!/usr/local/bin/julia


"""
  # Module NCload

  Load netCDF data with DSMACC output into dataframes for species concentrations
  and reaction rates.


  # Functions

  - get_ncdata (public)
  - ncload (private)
"""
module NCload

  ##################
  ###  PREAMBLE  ###
  ##################

  # Export public functions
  export get_ncdata

  # Assume either DSMACC main folder or DSMACC/AnalysisTools/DSMACCplot
  # as current directory, other wise add/adjust folder path here:
  if splitdir(pwd())[2] == "DSMACCplot"  def_dir = "../../save/results"
  else def_dir = "./save/results"
  end

  # Loading Julia and self-made modules
  # Define directories of self-made modules in main script
  # (absolute or relative paths to location, where main script is called)
  using PyCall, DataFrames
  using fhandle: test_file


###################
###  FUNCTIONS  ###
###################

"""
    get_ncdata(ncfile::String)

  Load netCDF data from `ncfile` and return dataframes with species concentrations
  `spec` and reaction rates `rate`.
"""
function get_ncdata(ncfile::String)

  # Check existance of file
  ncfile = test_file(ncfile,default_dir=def_dir)
  # Define path of python library and call it to read netCDF data
  unshift!(PyVector(pyimport("sys")["path"]), "py.lib")
  @pyimport ncdata

  # Read netCDF data from file
  nc = ncdata.get(ncfile)

  # Load data into dataframes for species and reaction rates
  spec = ncload(nc,"spec")
  rate = ncload(nc,"rate")

  # Return dataframes
  return spec, rate
end


"""
    ncload(nc::Dict{Any,Any},what::String="spec")

  Load data from dictonary `nc` with netCDF data from DSMACC output into a dataframe
  and return it. Only `what` is specified (species concentrations or reaction rates)
  will be returned (by default species concentrations).
"""
function ncload(nc::Dict{Any,Any},what::String="spec")

  # Initilise dataframe
  df = DataFrame()
  # Define array name with headers
  cols = what[1]*'c'
  # Loop over headers and add dataframe columns
  for head in nc[cols]
    df[Symbol(head)] = Float64[]
  end
  # Add data to dataframe
  for i = 1:length(nc[what][:,1])
    push!(df,nc[what][i,:])
  end

  # Return dataframe
  return df
end

end #module NCload
