#!/usr/local/bin/julia


"""
# Module ropa

Analyse sink and source fluxes from DSMACC output.
"""
module ropa


##################
###  PREAMBLE  ###
##################

println("initialise...")
# Define location of external self-made modules
# (Add or modify to include your own directories)
# Local Mac:
if isdir("/Applications/bin/data/jl.mod") &&
  all(LOAD_PATH.!="/Applications/bin/data/jl.mod")
  push!(LOAD_PATH,"/Applications/bin/data/jl.mod")
end
# earth0:
if isdir("~/Util/auxdata/jl.mod") && all(LOAD_PATH.!="~/Util/auxdata/jl.mod")
  push!(LOAD_PATH,"~/Util/auxdata/jl.mod")
end

# Add path of internal self-made modules
if isdir(joinpath(Base.source_dir(),"jl.mod")) &&
  all(LOAD_PATH.!=joinpath(Base.source_dir(),"jl.mod"))
  push!(LOAD_PATH,joinpath(Base.source_dir(),"jl.mod"))
end
if isdir(joinpath(Base.source_dir(),"jl.mod/ropa")) &&
  all(LOAD_PATH.!=joinpath(Base.source_dir(),"jl.mod/ropa"))
  push!(LOAD_PATH,joinpath(Base.source_dir(),"jl.mod/ropa"))
end

# Load modules/functions/python libraries
using PyPlot, DataFrames
using NCload, prepare_ropa

# Define the netCDF files from the first script argument
for i = 1:1-length(ARGS)  push!(ARGS,"")  end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

function get_fluxes(spc,rct_list,rates,fl)
  fluxes = Matrix{Any}(0,2)
  for i = 1:length(rct_list)
    try idx = find(rct_list[i][:,2].==spc)[1]
      if fl=="+"  fluxes = vcat(fluxes,[rates[i,1] rates[i,2]./rct[i][idx,1]])
      else        fluxes = vcat(fluxes,[rates[i,1] rates[i,2].*rct[i][idx,1]])
      end
    end
  end

  return fluxes
end


function prodloss(spc_arr,edct,prdct,flux_data)
  source_scen = []; sink_scen = []
  for s = 1:length(flux_data)
    sources = Matrix{Any}(length(spc_arr[s]),2)
    sinks = Matrix{Any}(length(spc_arr[s]),2)
    sources[:,1] = spc_arr[s]; sinks[:,1] = spc_arr[s]
    for (i, spc) in enumerate(spc_arr[s])
      println(spc)
      sources[i,2] = get_fluxes(spc,prdct[s],flux_data[s],"+")
      sinks[i,2] = get_fluxes(spc,edct[s],flux_data[s],"-")
    end
    push!(source_scen,sources); push!(sink_scen,sinks)
  end

  return source_scen, sink_scen
end


function flux_rates(ed_scen,specs,rates,rxn)

  flux_data = []
  for i = 1:length(ed_scen)
    flux_rate = Matrix{Any}(0,2)
    for j = 1:length(ed_scen[i])
      v = rates[i][Symbol(rxn[i][j])].*specs[i][Symbol(educt[i][j][1,2])].^educt[i][j][1,1]
      for k = 2:length(educt[i][j][:,2])
        v .*= specs[i][Symbol(educt[i][j][k,2])].^educt[i][j][k,1]
      end
      flux_rate = vcat(flux_rate, [rxn[i][j] [v]])
    end
    push!(flux_data,flux_rate)
  end

  return flux_data
end

function net_flux(educt,product,flux_data)

  for s = 1:length(educt)
    for i=1:length(educt[s])-1
      for j=i+1:length(product[s])
        if sort(educt[s][i][:,2])==sort(product[s][j][:,2]) && sort(product[s][i][:,2])==sort(educt[s][j][:,2])
          net = flux_data[s][i,2].-flux_data[s][j,2]
          flux_data[s][i,2] = zeros(length(net)); flux_data[s][j,2] = zeros(length(net))
          f=find(net.>0); b=find(net.<0)
          flux_data[s][i,2][f] = net[f]; flux_data[s][j,2][b] = net[b].⋅-1
          flux_data[s][i,1]=replace(flux_data[s][i,1],"-->","<-->>")
          flux_data[s][j,1]=replace(flux_data[s][j,1],"-->","<-->>")
        end
      end
    end
  end

  return rxn, flux_data
end


#####################
###  MAIN SCRIPT  ###
#####################

println("load data...")
# Load DSMACC output
specs = []; rates = []
for file in split(ARGS[1])
  spec, rate = NCload.get_ncdata(file)
  push!(specs,spec); push!(rates,rate)
end
(specs,rates)

# Get relevant species and reactions for each scenario as array of arrays
spc, rxn = get_names(specs,rates)
educt, product = split_rxn(rxn)
flux_data = flux_rates(educt,specs,rates,rxn)
flux_data = net_flux(educt,product,flux_data)
source, sink = prodloss(spc,educt,product,flux_data)

# Find sinks and sources for each species
# prodloss()


end #module ropa

#=
for i=1:length(educt[2])-1
  for j=2:length(educt[2])
    if sort(educt[2][i][:,2])==sort(product[2][j][:,2]) && sort(product[2][i][:,2])==sort(educt[2][j][:,2])
      println(rxn[2][i],'\n',rxn[2][j],"\n",educt[2][i],'\n',product[2][j],"\n")
    end
  end
end
=#
