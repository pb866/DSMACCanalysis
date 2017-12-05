__precompile__()


"""
# Module groupSPC

Module to group species and concentrations by the properties
- Chromophore class
- O:C ratio
- Molecule size

# Functions (public)

- readDB
- translateNMVOC
- translateSPC
- group_CC
- add_CC
"""
module groupSPC


  ##################
  ###  PREAMBLE  ###
  ##################

  # Export public functions
  export readDB,
         translateNMVOC,
         translateSPC,
         group_CC,
         add_CC

  # Loading external and internal self-made modules
  # Define directory of modules in main script
  # (absolute or relative paths to location, where main script is called)
  using DataFrames
  using fhandle


###################
###  FUNCTIONS  ###
###################

"""
    readDB(ifile::String)

From database file `ifile`, read different name definitions and molar masses of
all MCM species and dataframe `db` with the columns `mcm`, `smiles`, `inchi`,
`gecko`, and `mmol`.
"""
function readDB(ifile::String)

  # Read database file
  lines = rdinp(ifile,1)
  # Initilise dataframe columns
  mcm = String[]; smiles = String[]; inchi = String[]; gecko = String[]
  mmol = Float64[]
  # Loop over database and split lines into columns
  for line in lines
    push!(mcm,strip(split(line,"&")[1]))
    push!(smiles,strip(split(line,"&")[2]))
    push!(inchi,strip(split(line,"&")[3]))
    push!(gecko,strip(split(line,"&")[4]))
    push!(mmol,tryparse(Float64,split(line,"&")[5]).value)
  end

  # Save columns in database
  db = DataFrame()
  db[:mcm] = mcm; db[:smiles] = smiles; db[:inchi] = inchi; db[:gecko] = gecko
  db[:mmol] = mmol

  # Return database
  return db
end #function readDB


"""
    translateNMVOC(dframes,spcDB::DataFrames.DataFrame)

From array with dataframes of species concentrations and the translation database,
return an array of arrays with species names translated to GECKO-A nomenclature.
"""
function translateNMVOC(dframes,spcDB::DataFrames.DataFrame)

  # Initialise "outer" array (array of arrays)
  gspc = []
  # Loop over dataframes in array
  for df in dframes
    # Initilise "inner" array for GECKO-A names in each scenario
    gnames = String[]
    # Loop over dataframe headers with species names
    for name in string.(names(df))
      # Translate species from MCM to GECKO-A names
      idx = find(spc==name  for spc in spcDB[:mcm])
      if isempty(idx)
        println("translateNMVOC: $name not found. Species ignored for grouping.")
      elseif spcDB[:gecko][idx[1]][1] != 'G' && spcDB[:gecko][idx[1]] != "CH4"
        # Save names on successful find in inner array
        push!(gnames,spcDB[:gecko][idx[1]])
      end
    end
    # Save inner array in outer array
    push!(gspc,gnames)
  end

  # Return array with arrays of GECKO-A names
  return gspc
end #function translateNMVOC


"""
    translateSPC(specs,db,sym_in,sym_out)

Translate a list of species `specs` with the help of the translation database `db`
from the chosen language with keyword `sym_in` to the chosen language with keyword
`sym_out` and return list of translated species names `tspc`.
"""
function translateSPC(specs,db,sym_in,sym_out)

  # Initilise list of translated species names
  tspc = String[]
  # Loop over input species names
  for spec in specs
    # Find index of species in translation database
    idx = find(spc==spec  for spc in db[Symbol(sym_in)])[1]
    # Return and save species with this index in the chosen output language
    push!(tspc,string(db[Symbol(sym_out)][idx]))
  end

  # Return list of translated species names
  return tspc
end #function translateSPC


"""
    group_CC(gspc,spcDB::DataFrames.DataFrame)

From an array with list of GECKO-A names `gspc` and the translation database `spcDB`,
group species into chromophore classes and return a list `CC` with the chromophore
class names and an array of arrays `chrom_class` with the GECKO-A names grouped by
chromophore classes for each scenario.

The following chromophore classes exist:
- `Ald`: Mono-aldehydes
- `Ket`: Mono-ketones
- `DiCar`: Di- and polycarbonyls
- `Kete`: Ketenes
- `Nitro`: Mono-nitrocompounds
- `DiNitro`: Di- and poly-nitrocompounds
- `Nit`: Mono-alkyl nitrates
- `DiNit`: Di- and poly-alkyl nitrates
- `PNit`: Mono-alkyl pernitrates
- `DiPNit`: di- and poly-alkyl pernitrates
- `PAN`: peroxyacyl nitrates
- `ROOH`: Mono-alkyl hydroperoxides
- `DiROOH`: Di- and poly-alkyl hydroperoxides
- `PAA`: Organic peroxy acids
- `SCI`: Stabilised Criegee intermediates
- `Poly`: Compounds with polyfunctional chromophores
- `NoChr`: Stable organic compounds with no chromophores
- `Rad`: organic radical compounds (excluding SCI)
- `Inorg`: Inorganic compounds
- `New`: New generic model species
"""
function group_CC(gspc,spcDB::DataFrames.DataFrame)

  # Initilise "outer" array for scenarios
  chrom_class = []
  # Loop over scenarios
  for data in gspc
    # Initilise arrays for different chromophore classes
    Ald = String[]; Ket = String[]; DiCar = String[]; Kete = String[];
    Nitro = String[]; DiNitro = String[]; Nit = String[]; DiNit = String[];
    PNit = String[]; DiPNit = String[]; PAN = String[];
    ROOH = String[]; DiROOH = String[]; PAA = String[]; SCI = String[];
    Poly = String[]; NoChr = String[]; Rad = String[]; New = String[]

    # Loop over species in each scenario
    for spc in data
      # Find all chromophores in chemical formula
      mspc  = translateSPC([spc],spcDB,"gecko","mcm")[1]
      ald   = length(matchall(r"(?<!-)CHO(?!(\(|-))",spc))
      ket   = length(matchall(r"(?<!-)C[0-9]*O(?!(\(|-))",spc))
      kete  = length(matchall(r"CdO",spc))
      nitro = length(matchall(r"\(NO2\)",spc))
      nit   = length(matchall(r"\(ONO2\)",spc))
      pnit  = length(matchall(r"(?<!CO)\(OONO2\)",spc))
      pan   = length(matchall(r"CO\(OONO2\)",spc))
      rooh  = length(matchall(r"(?<!CO)\(OOH\)",spc))
      paa   = length(matchall(r"CO\(OOH\)",spc))
      sci   = length(matchall(r"\.\(OO\.\)(?!\*)",spc))
      # Save number of chromophores in each class in an array
      chrom = [ald,ket,kete,nitro,nit,pnit,pan,rooh,paa,sci]
      # Double check correctness of GECKO-A string
      if uppercase(spc[1])!='C' && spc[1]!='-'  &&  spc[1:3]!="new"
        println("Warning! Species $spc does not start with carbon group, ether group")
        println("or key 'new' for new generic model species. Unpredicted results possible.")
      end

      # Save MCM names in the respective chromophore classes
      if     spc[1:3] == "new"                                    push!(New,mspc)
      elseif length(matchall(r"\.",spc))>0 &&
             length(matchall(r"\.\(OO\.\)(?!\*)",spc))==0         push!(Rad,mspc)
      elseif all(chrom.==0)                                       push!(NoChr,mspc)
      elseif all(chrom[2:end].==0) && chrom[1]==1                 push!(Ald,mspc)
      elseif chrom[1]==0 && all(chrom[3:end].==0) && chrom[2]==1  push!(Ket,mspc)
      elseif (chrom[1]>0 || chrom[2]>0) && all(chrom[3:end].==0)  push!(DiCar,mspc)
      elseif all(chrom[1:2].==0) && all(chrom[4:end].==0)         push!(Kete,mspc)
      elseif all(chrom[1:3].==0) && all(chrom[5:end].==0)
        if chrom[4]==1                                            push!(Nitro,mspc)
        else                                                      push!(DiNitro,mspc)
        end
      elseif all(chrom[1:4].==0) && all(chrom[6:end].==0)
        if chrom[5]==1                                            push!(Nit,mspc)
        else                                                      push!(DiNit,mspc)
        end
      elseif all(chrom[1:5].==0) && all(chrom[7:end].==0)
        if chrom[6]==1                                            push!(PNit,mspc)
        else                                                      push!(DiPNit,mspc)
        end
      elseif all(chrom[1:6].==0) && all(chrom[8:end].==0)         push!(PAN,mspc)
      elseif all(chrom[1:7].==0) && all(chrom[9:end].==0)
        if chrom[8]==1                                            push!(ROOH,mspc)
        else                                                      push!(DiROOH,mspc)
        end
      elseif all(chrom[1:8].==0) && all(chrom[9:end].==0)         push!(PAA,mspc)
      elseif all(chrom[1:9].==0)                                  push!(SCI,mspc)
      else                                                        push!(Poly,mspc)
      end
    end
    # Save an array with the chromophore classes for each scenario
    push!(chrom_class,[Ald,Ket,DiCar,Kete,Nitro,DiNitro,Nit,DiNit,PNit,DiPNit,PAN,
      ROOH,DiROOH,PAA,SCI,Poly,NoChr,Rad,New])
  end
  # Save a list with the chromophore class names
  CC = ["Ald","Ket","DiCar","Kete","Nitro","DiNitro","Nit","DiNit","PNit","DiPNit",
        "PAN","ROOH","DiROOH","PAA","SCI","Poly","NoChr","Rad","New"]

  # Return the list of class names and
  # the array with chromophore classes for each scenario
  return CC, chrom_class
end #function group_CC


"""
    add_CC(specs,CC,classes)

Append the array with dataframes of species concentrations for each scenario `specs`
by the concentrations of all species of the various chromophore `classes`
(with the names defined in the list `CC`).
"""
function add_CC(specs,CC,classes)

  # Loop over scenarios
  for i = 1:length(CC)
    # Loop over chromophore classes
    for (j, class) in enumerate(classes)
      # Initilise concentration arrays for each chromophore class
      specs[i][Symbol(class)] = zeros(length(specs[i][1]))
      # Loop over species in each chromophore class
      for spc in CC[i][j]
        # Add concentrations of individual species to chromophore class
        specs[i][Symbol(class)] += specs[i][Symbol(spc)]
      end
    end
  end

  # Return appended array
  return specs
end #function add_CC

end #module groupSPC
