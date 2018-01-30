__precompile__()


"""
# Module jlplot

Define type of plots and species/reactions to be plotted from an input file.
Generate arrays with names of DSMACC netCDF files, scenario names,
species/reactions to be plotted and further indices needed for the processing
of the data.


# Functions

## Public
- commission_plot
- get_scenario
- get_settings
- prepare_plots
- get_stackdata

## Private
- rm_blanklines
"""
module jlplot


##################
###  PREAMBLE  ###
##################

# Export public functions
export commission_plot,
       get_scenario,
       get_settings,
       prepare_plots,
       DSMACCoutput,
       get_stackdata,
       def_night

# Find directory of module source code
cdir=Base.source_dir()

# Loading external and internal self-made modules
# Define directory of modules in main script
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
# Current directory
if all(LOAD_PATH.!=cdir)  push!(LOAD_PATH,cdir)  end

using fhandle
using make_plots: night_shade
using NCload
using DataFrames


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    commission_plot(ifile::String="plot.inp")

Read input file `ifile` and return Array with lines `commission` as well as
indices for the beginning of the scenario section `scen`, the beginning of the
settings sections `sett`, and the beginning and end of the plotting section
`beg_plt` and `end_plt`.
"""
function commission_plot(ifile::String="plot.inp")
  # Assume either DSMACC main folder or DSMACC/AnalysisTools/DSMACCanalysis
  # as current directory, other wise add/adjust folder path here:
  if splitdir(pwd())[2] == "DSMACCanalysis"  def_dir = "."
  else def_dir = "./AnalysisTools/DSMACCanalysis"
  end
  # Read and store lines from input file
  commission = rdinp(ifile, default_dir=def_dir)
  commission = strip.(commission)
  # Remove comments in file content
  for i = length(commission):-1:1
    if startswith(commission[i],"#")  deleteat!(commission,i)
    else try
      commission[i] = replace(commission[i],match(r"#.*",commission[i]).match,"")
    end end
  end

  # Find indices for beginning of scenario definitions
  # and beginning/end of plotting section
  scen = find([contains(line,"Scenarios:") for line in commission])[1] + 1
  sett = find([contains(line,"Settings:") for line in commission])
  beg_plt = find([contains(line,"Plotting:") for line in commission])[1] + 1
  end_plt = find([contains(line,"Comments:") for line in commission])
  # Correct end of plotting section if not followed by Comments section
  if isempty(sett)  sett = 0
  else              sett = sett[1] + 1
  end
  # Correct end of plotting section if not followed by Comments section
  if isempty(end_plt)
    end_plt = length(commission)
  else
    end_plt = end_plt[1] - 1
    commission = commission[1:end_plt]
  end
  # Compile indices in an array
  fidx = [scen, sett, beg_plt, end_plt]

  # Remove blank lines (whitespaces allowed) between section headings
  # and first content line
  commission, fidx = rm_blanklines(commission,fidx)

  # Return file content and indices
  return commission, fidx
end #function commission_plot


"""
    get_scenario(lines,strt)

From `lines` in input file and index `strt` for beginning of the Scenario section,
retrieve and return arrays `ncfile` and `label` with the directory and name of the
DSMACC netCDF files and names for the respective scenarios.

If no scenario names are specified, the netCDF file names without the extension
will be used.
"""
function get_scenario(lines,strt)

  # Get nc file names
  lines[strt] = replace(lines[strt],r",|;"," ") # allow other separators than spaces
  ncfile = split(lines[strt])

  # Get scenario names
  label = String[]
  if lines[strt+1]==""
    # Get scenario names from file names, if no labels are given in input file
    for scen in ncfile  push!(label,basename(splitext(scen)[1]))  end
  else
    # Read scenario names from file (labels must be wrapped in double quotes)
    s = Int64[]; si = 1
    for i = 1:length(ncfile)*2
      push!(s,searchindex(lines[strt+1],"\"",si))
      si = s[end] + 1
    end
    for i = 1:2:length(ncfile)*2
      push!(label,lines[strt+1][s[i]+1:s[i+1]-1])
    end
  end

  # Error handling
  if length(label) != length(ncfile)
    # Stop if number of labels is unequal to number of nc files
    println("There must be exactly as much labels define in the input file as nc file")
    println("or no labels at all. (Don't use double quotes in labels.) Script stopped.")
    exit()
  end
  for lab in label
    # Stop if duplicate labels are defined
    len = length(find(s==lab for s in label))
    if len > 1  println("Error: Non-unique labels! Script stopped."); exit()  end
  end

  # Return arrays with file names and labels
  return ncfile, label
end #function get_scenario


"""
    get_settings(lines,sett_idx)

From the `lines` in the input plot file and the index for the start of the
settings section in the lines `sett_idx`, find and return the lower and upper
cut-off parameters for minor/major fluxes (as array `[llim, ulim]`), the switch
to calculate (or not) net cycles of chemical fluxes, and the colour and
transparency of night-time shading in plots or return the default values of
llim = 0.05, ulim = 0.7, cycles = "reduce", and no night-time shading.
"""
function get_settings(lines,sett_idx)

  # Set default cut-off parameters
  mcm = "v3.3.1"
  llim = 0.05
  ulim = 0.7
  cycles = "reduce"
  nightcol = "w"; ntrans = 0.0
  t_frmt = "TIME"
  fig = "off"
  # Overwrite parameters with values from the Settings section, if defined
  if sett_idx!=0
    i = sett_idx
    while lines[i] != ""
      if lines[i][1:4] == "MCM:"
        mcm = strip(lines[i][5:end])
      elseif lines[i][1:5]=="time:"
        t_frmt = strip(lines[i][6:end])
      elseif lines[i][1:6]=="night:"
        nightcol, ntrans = strip.(split(lines[i][7:end],"/"))
        ntrans = float(ntrans)
      elseif lines[i][1:8]=="cut-off:"
        lines[i] = replace(lines[i],r",|;"," ") # allow other separators than spaces
        llim, ulim = float.(split(lines[i][9:end]))
      elseif lines[i][1:7]=="cycles:"
        cycles = strip(lines[i][8:end])
      elseif lines[i][1:4]=="Fig:"
        fig = strip(lines[i][6:end])
      end
      i += 1
    end
  end

  # Return lower and upper cut-off
  return [llim, ulim], cycles, [nightcol, ntrans], mcm, t_frmt, fig
end #function get_settings


"""
    prepare_plots(commission,label)

From the plotting section of the input file `commission` and the `label`s of
each scenario, retrieve the index `icase` of the scenario to plot, `what` to
plot (species concentrations/reaction rates/fluxes) and the `plotdata` with the
species/reactions that go into each plot.
"""
function prepare_plots(commission,label)

  # Initilise
  new_case = true
  icase = []
  what = []
  unit = []
  plotdata = []
  pdata = []
  # Loop over plotting section of input file
  for (i, line) in enumerate(commission)
    if line == "" && commission[i-1]==""
      continue
    elseif line == ""
      # On empty lines, reset flag for new case, collect current plot data in array
      # and reset temporary memory for plot data
      new_case = true
      push!(plotdata,pdata)
      pdata = []
    elseif new_case
      # First line or lines after empty lines look for plot type definition
      new_case = false #set flag for new case to false
      # split line into pre- and post-colon
      icol = searchindex(line,":")
      # Get post-colon data and save what to plot (concentrations/rates) and units
      wh = ""; un = ""
      try wh, un = split(line[icol+1:end],"/")
      catch
        # If units are obsolete, use mlc·cm-3 (s-1) as default
        wh = strip(line[icol+1:end])
        un = "mlc"
      end
      push!(what,strip(wh)); push!(unit,strip(un))

      # Look for the scenarios to plot before the colon
      # Read scenario labels
      idx = Int64[]
      line = replace(line,r",|;"," ") # allow other separators than spaces
      lab = split(line[1:icol-1])
      # Get indices of each scenario
      for l in lab
        ic = find(s==l for s in label)[1]
        push!(idx, ic)
      end
      # save indices in array icase
      push!(icase,idx)
    else
      # Save array with species/reactions to be plotted to array plotdata
      # when reading general data lines
      # Allow commas, semicolons or whitespace as separator
      line = replace(line,r",|;"," ")
      push!(pdata,strip.(split(line)))
    end
  end

  if commission[end] != ""  push!(plotdata,pdata)  end

  # Return final arrays
  return icase, what, unit, plotdata
end #function prepare_plots


"""
    get_stackdata(spc_list,case,ydata,unit)

Plot the concentrations of the species in `spc_list` for the current `case` using
the `ydata` in the specified `unit` and return y data in the correct unit for the
boundaries and the areas in the graph.
"""
function get_stackdata(spc_list,case,ydata,unit)
  # Initialise output arrays
  ystack = DataArray[]; ylines = DataArray[]
  # Define areas for stack plot
  for spc in spc_list
    push!(ystack,ydata[case][Symbol(spc)])
  end
  # Define boundary lines for stack plot
  for i = 1:length(ystack)
    push!(ylines,sum(ystack[1:i]))
  end

  # Perform unit conversions
  if unit=="ppm"
    for l = 1:length(ylines)
      ylines[l] .*= 1.e6./ydata[case][:M]
      ystack[l] .*= 1.e6./ydata[case][:M]
    end
  elseif unit=="ppb"
    for l = 1:length(ylines)
      ylines[l] .*= 1.e9./ydata[case][:M]
      ystack[l] .*= 1.e9./ydata[case][:M]
    end
  elseif unit=="ppt"
    for l = 1:length(ylines)
      ylines[l] .*= 1.e12./ydata[case][:M]
      ystack[l] .*= 1.e12./ydata[case][:M]
    end
  end

  # Return boundaries and areas to be plotted over time
  return ylines, ystack
end


"""
    def_night(rates, ntrans)

From the array `rates` with DataFrames of all scenarios holding the reaction rates,
return a n×2 matrix indices for the model times at the beginning/end of night for
all scenarios. If the transparency value for night-time shading `ntrans` is set
to 0 in the settings, skip the routine and return an empty array.
"""
function def_night(rates, ntrans)
  # Initilise array
  night = Any[]
  # Skip, if night-time shading is switched off in settings
  if ntrans != 0.0
    # Loop over scenarios
    for (i, scen) in enumerate(rates)
      # Get indices for times of sunsets/sunrises
      nght = Matrix{Int64}(0,2)
      try nght = night_shade(scen[Symbol("NO2-->O+NO")])
      catch
        try nght = night_shade(scen[Symbol("NO2-->NO+O")])
        catch
          print("J(NO2) not found in rates. ")
          println("Night-time shading omitted in Scenario $i.")
        end
      end
      # Save starts/ends of all nights for each scenario
      push!(night, nght)
    end
  end
  if isempty(night)  night = ones(Float64,length(rates)).⋅(-1)  end

  # Return array with bounderies of all nights in all scenarios
  return night
end


###########################
###  PRIVATE FUNCTIONS  ###
###########################

"""
    rm_blanklines(lines,idx)

Starting at index `idx` in the string array `lines` remove all blank lines
until the first line with non-whitespace characters is found.
"""
function rm_blanklines(lines,idx)
  # Loop over different sections
  for i = length(idx)-1:-1:1
    # Skip for missing 2nd section
    if idx[i] == 0  continue  end
    # Delete leading blank lines until first line with content is found
    while strip(lines[idx[i]])==""
      deleteat!(lines,idx[i])
      idx[i+1:end] -= 1
    end
  end
  # Revert index correction for missing second section
  idx[2] = max(0,idx[2])

  # Return adjusted array with lines
  return lines, idx
end #function rm_blanklines

end #module jlplot
