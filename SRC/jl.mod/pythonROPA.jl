

module pythonROPA

# Define location of external self-made modules
# (Add or modify to include your own directories)
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

# Add path of internal self-made modules
ropa_dir = joinpath(Base.source_dir(),"ropa")
if all(LOAD_PATH.!=ropa_dir)  push!(LOAD_PATH,ropa_dir)  end

# Load Julia modules
using DataFrames, PyCall
using plot_ropa
# Load python libraries
unshift!(PyVector(pyimport("sys")["path"]),
  normpath(joinpath(Base.source_dir(),"../py.lib")))
@pyimport matplotlib.backends.backend_pdf as pdf
@pyimport ropa_fn


### USER SPECIFICATIONS ###
# Set DSMACC output time step in seconds
dt = 600
spc = [["OH"],["HO2"],["O3"],["NO"],["NO2"],["NO","NO2","N2O5"]]
fct = [[1],[1],[1],[1],[1],[1,1,2]]
scen = ["DATA/TOLm3.nc","DATA/TOLm4.nc"]

# Initialise sources and sinks array
sources = Vector{Any}(length(scen)); sinks = Vector{Any}(length(scen))
for i = 1:length(scen)
  sources[i] = Matrix{Any}(0,2); sinks[i] = Matrix{Any}(0,2)
end

# Add selected species to sources and sinks matrix
for i = 1:length(spc) #loop over species
  # Call python ROPA version
  ropa_src, ropa_snk, src_names, snk_names = ropa_fn.prodloss(scen,spc[i],fct[i])
  # ropa plotting routines assume positive sink fluxes, converte them to positive floats
  ropa_snk .*= -1

  for s = 1:length(ropa_src)
    # Save data to arrays of DataFrames with an DataFrame for each scenario
    # Separate arrays for sources and sinks
    src = Matrix{Any}(0,2); snk = Matrix{Any}(0,2)
    for rxn = 1:length(ropa_src[s][:,1])
      src = vcat(src, [String(src_names[s][rxn]) [ropa_src[s][rxn,:]]])
    end
    for rxn = 1:length(ropa_snk[s][:,1])
      snk = vcat(snk, [String(snk_names[s][rxn]) [ropa_snk[s][rxn,:]]])
    end

    sources[s] = vcat(sources[s],[join(spc[i],"+") [src]])
    sinks[s] = vcat(sinks[s],[join(spc[i],"+") [snk]])
    #=
    srcdf = DataFrame(); snkdf = DataFrame()
    srcdf[:time] = collect(1/3600*dt:1/3600*dt:length(src[s][1,:])/3600*dt)
    snkdf[:time] = collect(1/3600*dt:1/3600*dt:length(snk[s][1,:])/3600*dt)
    for (i, name) in enumerate(Symbol.(String.(src_names[s])))
      srcdf[name] = src[s][i,:]
    end
    for (i, name) in enumerate(Symbol.(String.(snk_names[s])))
      snkdf[name] = snk[s][i,:]
    end
    push!(sources,srcdf)
    push!(sinks,snkdf)
    =#
  end

end #loop over species

# Calculate model times for each scenario from the specified time step
tdf = Vector(length(scen))
t = DataFrame(TIME = collect(dt:dt:length(sources[1][1,2][1,2])⋅dt))
tdf .= t

# Get single species for plotting
spec = String[]
for s in spc
  if length(s)==1  push!(spec,s[1])  end
end
# Set scenario titels if not defined in user specifications
if !isdefined(:titel)
  titel = []
  for fname in scen
    push!(titel,Base.splitext(basename(fname))[1])
  end
end
# Plot data
pdffile = pdf.PdfPages("PETE.pdf")
for species in spec  for s = 1:length(scen)
  modtime, Ysrc, Ysnk, Ysrc_rev, Ysnk_rev =
    load_plotdata(species,s,sources,sinks,tdf,titel)# Output flux plots

  fig = plot_data(species,titel[s],modtime,Ysrc,Ysnk)
  if fig != nothing  pdffile[:savefig](fig)  end
  # if fig != nothing  fig[:show]()  end
  # Output revised flux plots, if major fluxes have been removed
  fig = plot_data(species,titel[s],modtime,Ysrc_rev,Ysnk_rev)
  if fig != nothing  pdffile[:savefig](fig)  end
end  end
# Close multipage pdf file
pdffile[:close]()

end #module pythonROPA
