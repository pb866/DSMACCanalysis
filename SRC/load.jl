### Read input plot file and netCDF files and define plot commission ###

println("load data...")
# Read input file
commission, fidx = commission_plot(ARGS[1])
# Get nc file names and scenario names
ncfiles, label = get_scenario(commission,fidx[1])
# Set cut-off for flux plots, switch to calculate inorganic net cycles,
# switch for night-time shading in plots and MCM version used in DSMACC
lims, cycles, pltnight, mcm, t_frmt, sfig = get_settings(commission,fidx[2])
# Read from input file, which plots to generate
icase, what, unit, plotdata = prepare_plots(commission[fidx[3]:fidx[4]], label)
# Save DSMACC output using a python script to read nc files and save in Julia format
specs, rates = DSMACCoutput(ncfiles)
# Data is stored in a dictionary distinguishing between "specs" and "rates"
# In each dictionary entry is is an array for each scenario (nc file)
# with dataframes for the respectives species concentrations or reaction rates

# Add sza to specs and rates dataframes
# specs, rates = addSZA(specs, rates)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

# Add folder for single pdfs
if sfig=="on"
  # Confirm overwrite of old data
  if ispath(normpath(joinpath(Base.source_dir(),"../FIG")))
    println("Overwrite plots in folder FIG?")
    println("(\'yes\': overwrite / \'<ENTER>\': quit)")
    confirm = input()
    if lowercase(confirm[1]) == 'y'
      rm("FIG",recursive=true)
    else
      println("exit code 99.")
      exit(99)
    end
  end
  # Create folder for single pdfs
  mkdir("FIG")
end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

### ROPA analysis ###

if any(what.=="fluxes")
  # Perform ROPA analysis for flux plots
  sources, sinks, concs = ropa(cycles=cycles, specs=specs, rates=rates)
  # Combine species for flux plots in single vector
  for i = 1:length(what)
    if what[i] == "fluxes"
      pldata = String[]
      for j = 1:length(plotdata[i])  for spc in plotdata[i][j]
        push!(pldata,spc)
      end  end
      plotdata[i] = pldata
    end
  end
end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

### Group species by properties

println("analyse data...")
# Load database with different species names
dbfile = normpath(joinpath(Base.source_dir(),"../DATA/MCMv331species_corrected.db"))
spcDB = readDB(dbfile)
# Translate species names from MCM to GECKO-A
gspc = translateNMVOC(specs,spcDB,MCM=mcm)
# Group species by properties
CC, OC, CN, chrom_class, OCratio_class, size_class = group_specs(gspc,spcDB,MCM=mcm)
# Add new classes to dataframes
specs = add_conc(specs,chrom_class,OCratio_class,size_class,CC,OC,CN)
