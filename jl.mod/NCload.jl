#!/usr/local/bin/julia


"""
# Module NCload

Load netCDF data with DSMACC output into dataframes for species concentrations
and reaction rates.


# Functions

## Public
- DSMACCoutput
- get_ncdata
- addSZA

## Private
- ncload
"""
module NCload

##################
###  PREAMBLE  ###
##################

# Export public functions
export DSMACCoutput,
       get_ncdata,
       addSZA

# Assume either DSMACC main folder or DSMACC/AnalysisTools/DSMACCanalysis
# as current directory, other wise add/adjust folder path here:
if splitdir(pwd())[2] == "DSMACCplot"  def_dir = "../../save/results"
else def_dir = "./save/results"
end

# Loading Julia and self-made modules
# Define directories of self-made modules in main script
# (absolute or relative paths to location, where main script is called)
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
using PyCall, DataFrames
using fhandle: test_file


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    DSMACCoutput(ncfiles)

Compile data of netCDF files in Julia readable formats.

Data will be stored in the arrays `specs` and `rates` holding dataframes with
species concentrations and reaction rates from each scenario in the `ncfiles`.
"""
function DSMACCoutput(ncfiles)

  # Assume either DSMACC main folder or DSMACC/AnalysisTools/DSMACCanalysis
  # as current directory, other wise add/adjust folder path here:
  if splitdir(pwd())[2] == "DSMACCanalysis"  def_dir = "../../save/results"
  else def_dir = "./save/results"
  end

  # Loop over all nc files, check for existance and read content
  # separately for species concentrations and reaction rates
  specs = []; rates = []
  for ncfile in ncfiles
    # Assign default path, if no path is specified for the nc files
    if dirname(ncfile)==""  ncfile = normpath(joinpath(def_dir,ncfile))  end
    # Read species concentrations/reaction rates of current nc file
    spc, rat = get_ncdata(ncfile)
    # Save concentrations/rates in array with all scenarios
    push!(specs,spc); push!(rates,rat)
  end

  # Return concentrations/rates of all scenarios
  return specs, rates
end #function DSMACCoutput


"""
    get_ncdata(ncfile::String)

Load netCDF data from `ncfile` and return dataframes with species concentrations
`specs` and reaction `rates`.
"""
function get_ncdata(ncfile::AbstractString)

  # Check existance of file
  ncfile = test_file(ncfile,default_dir=def_dir)
  # Define path of python library and call it to read netCDF data
  unshift!(PyVector(pyimport("sys")["path"]),
    joinpath(Base.source_dir(),"py.lib"))
  @pyimport ncdata

  # Read netCDF data from file
  specs, rates = ncdata.load(ncfile)

  # Transform Python Pandas DataFrames to Julia DataFrames
  specs, rates = ncload(specs, rates)

  # Return dataframes
  return specs, rates
end


"""
    addSZA(specs, rates)

Add columns `sza` and `SZA` with solar zenith angles in deg and rad, respecitively,
to arrays `specs` with species concentrations and array `rates` with reaction
rates for each scenario.
"""
function addSZA(specs, rates)

  # Loop over scenarios
  for s = 1:length(specs)
    # Get Julian day and time of the day in seconds (jsec)
    jday = Dates.dayofyear.(specs[s][:JTIME])
    # jd = Dates.value.(Dates.Second.(Dates.Day.(jday)))
    jsec = Dates.value.(Dates.Second.(Dates.Hour.(specs[s][:JTIME]))) .+
           Dates.value.(Dates.Second.(Dates.Minute.(specs[s][:JTIME]))) .+
           Dates.value.(Dates.Second.(specs[s][:JTIME]))

    # Get latitude in radian from model output
    lat = deg2rad.(specs[s][:LAT])

    # Calculate hourly angles in radians
    w = deg2rad.(15.0*(jsec./3600-12.))

    # Calculate declination
    φ = 2.0π.⋅(jday-1)/365
    dec = (0.006918 - 0.399912cos.(φ) + 0.070257sin.(φ)
          - 0.006758cos.(2.0φ) + 0.000907sin.(2.0φ)
          - 0.002697cos.(3.0φ) + 0.00148sin.(3.0φ))

    # Calculate sza
    χ = acos.(sin.(dec) .⋅ sin.(lat) .+ cos.(dec) .⋅ cos.(lat) .⋅ cos.(w))
    sza = rad2deg.(χ)
    # Save sza columns in rad and deg in specs and rates dataframes
    specs[s][:SZA] = χ
    specs[s][:sza] = sza
    rates[s][:SZA] = χ
    rates[s][:sza] = sza
  end

  # Return ammended dataframe arrays
  return specs, rates
end #function calcSZA


###########################
###  PRIVATE FUNCTIONS  ###
###########################

"""
    ncload(spc::PyCall.PyObject,rat::PyCall.PyObject)

Load data from Python objects `spc` and `rat` with species concentrations and
reaction rates from the DSMACC netCDF model output into a dataframe
and return it. Correct the first time step of the model time in column `JTIME`
and add a column `TIME` with model times starting at `0`.
"""
function ncload(spc::PyCall.PyObject,rat::PyCall.PyObject)

  # Define column header for Julia specs and rates DataFrame
  # from python pandas dataframes
  spc_cols = Symbol[]
  for s in spc[:columns]  push!(spc_cols, Symbol(string(s)))  end
  rat_cols = Symbol[]
  for s in rat[:columns]  push!(rat_cols, Symbol(string(s)))  end

  # Get Julian and model time (starting at 0)
  jtime = DateTime[]
  for t in spc[:TIME][:index]  push!(jtime,t)  end
  dt = Dates.Period(jtime[2]-jtime[1])
  dt=Dates.Second(dt)
  dt = convert(Float64,Dates.value(dt))/3600
  t = 0.
  modtime = Float64[]
  for i = 1:length(jtime)  push!(modtime,t); t += dt  end

  # Initilise Julia DataFrames with times
  specs = DataFrame(JTIME = jtime, TIME = modtime)
  rates = DataFrame(JTIME = jtime, TIME = modtime)
  # Add data to DataFrames
  for s in spc_cols[2:end]  specs[s] = spc[s][:values]  end
  for s in rat_cols[2:end]  rates[s] = rat[s][:values]  end


  # Return dataframe
  return specs, rates
end

end #module NCload
