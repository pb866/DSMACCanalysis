__precompile__()


"""
# Module groupSPC

Module to group species and concentrations by the properties
- Chromophore class
- O:C ratio
- Molecule size

# Functions

## public
- readDB
- translateNMVOC
- translateSPC
- group_specs
- add_conc

## private
- init_chrom
- init_ChromClass
- group_CC
- init_OCclass
- group_OC
- init_SizeClass
- group_CN
- group_conc
"""
module groupSPC


  ##################
  ###  PREAMBLE  ###
  ##################

  # Export public functions
  export readDB,
         translateNMVOC,
         translateSPC,
         group_specs,
         add_conc

  # Loading external and internal self-made modules
  # Define directory of modules in main script
  # (absolute or relative paths to location, where main script is called)
  using DataFrames
  using fhandle


##########################
###  PUBLIC FUNCTIONS  ###
##########################

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
    group_specs(gspc,spcDB::DataFrames.DataFrame)

From an array with list of GECKO-A names `gspc` and the translation database `spcDB`,
group species into classes of different chromophore types, O:C ratio ranges, and
molecule sizes, respectively.

Return a list `CC` with the respective class names for each property and arrays
of arrays `chrom_class`, `OCratio_class`, and `size_class` with the GECKO-A names
grouped by the respective properties for each scenario.


### Chromophore classes
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

### O:C ratio classes
- `oc1`:       O:C ≤ 0.5
- `oc2`: 0.5 < O:C ≤ 1.0
- `oc3`: 1.0 < O:C ≤ 2.0
- `oc4`: 2.0 < O:C

### Size classes
- `sc1`: CN = 1
- `sc2`: CN = 2
- `sc3`: CN = 3
- `sc4`: CN = 4
- `sc5`: CN = 5
- `sc6`: CN = 6
- `sc7`: CN = 7
- `sc8`: CN = 8
- `sc9`: CN = 9
- `sc0`: CN ≥ 10

CN stands for carbon number, although ether groups count towards the molecule's
size as well.

"""
function group_specs(gspc,spcDB::DataFrames.DataFrame)

  # Initilise "outer" array for scenarios
  chrom_class = []; OCratio_class = []; size_class = []
  # Loop over scenarios
  for data in gspc
    # Initilise arrays for different chromophore classes
    chromophores = init_ChromClass()

    # Initilise arrays for different O:C ratio classes
    OCclasses = init_OCclass()
    # Initilise arrays for different size classes (0: size ≥ 10)
    SizeClasses = init_SizeClass()

    #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#


    # Loop over species in each scenario
    for spc in data
      # Double check correctness of GECKO-A string
      if uppercase(spc[1])!='C' && spc[1]!='-'  &&  spc[1:3]!="new"
        println("Warning! Species $spc does not start with carbon group, ether group")
        println("or key 'new' for new generic model species. Unpredicted results possible.")
      end
      ### Initilise ###
      # Translate species name back to MCM nomenclature
      mspc  = translateSPC([spc],spcDB,"gecko","mcm")[1]

      # Determine all chromophores in the molecule
      chrom = init_chrom(spc)
      # Determine molecule's size (number of carbon and ether groups)
      CN=length(matchall(r"C(?!l)"i,spc))+length(matchall(r"-O"i,spc))
      # Determine O:C ratio
      # O: number of explicitly written 'O' + implied O using abreviations as in NO2
      O=length(matchall(r"O"i,spc))+length(matchall(r"O2\)"i,spc))
      # C: number of C, but no C in Cl + number of aromatic c (therefore case-insensitive)
      C=length(matchall(r"C(?!l)"i,spc))
      # Final ratio:
      OCratio=O/C

      #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

      # Save MCM names in the respective chromophore classes
      chromophores = group_CC(spc, mspc, chrom, chromophores)

      # Save MCM names in the respective O:C ratio classes
      OCclasses = group_OC(OCclasses,OCratio,spc,mspc)

      # Save MCM names in the respecitive size classes
      SizeClasses = group_CN(SizeClasses,CN,spc,mspc)
    end
    # Save an array with the chromophore, O:C ratio, and size classes for each scenario
    push!(chrom_class,chromophores)
    push!(OCratio_class,OCclasses)
    push!(size_class,SizeClasses)
  end
  # Save a list with the chromophore class names
  CC = ["Ald","Ket","DiCar","Kete","Nitro","DiNitro","Nit","DiNit","PNit","DiPNit",
        "PAN","ROOH","DiROOH","PAA","SCI","Poly","NoChr","Rad","New"]
  OC = ["oc1","oc2","oc3","oc4"]
  CN = ["sc1","sc2","sc3","sc4","sc5","sc6","sc7","sc8","sc9","sc0"]

  # Return the list of class names and
  # the array with chromophore classes for each scenario
  return CC, OC, CN, chrom_class, OCratio_class, size_class
end #function group_specs


"""
    add_conc(specs,chrom_class,OCratio_class,size_class,CC,OC,CN)

Append the array with dataframes of species concentrations for each scenario in
`specs` by the concentrations of all species with certain properties.

Properties include classes by
- chromophore type (`chrom_class`)
- O:C ratio (`OCratio_class`)
- molecule's size (`size_class`)

Class names for the dataframe headers are defined in the arrays `CC`, `OC`, and
`CN`, respecitively.
"""
function add_conc(specs,chrom_class,OCratio_class,size_class,CC,OC,CN)

  # Loop over scenarios
  for i = 1:length(chrom_class)
    # Add chromophore classes to DSMACC dataframes
    specs[i] = group_conc(specs[i],chrom_class[i],CC)
    # Add O:C reactions classes to DSMACC dataframes
    specs[i] = group_conc(specs[i],OCratio_class[i],OC)
    # Add size classes to DSMACC dataframes
    specs[i] = group_conc(specs[i],size_class[i],CN)
  end

  # Return appended array
  return specs
end #function add_conc (add_conc)


###########################
###  PRIVATE FUNCTIONS  ###
###########################

"""
    init_chrom(spc)

Return an array `chrom` with the number of each chromophore type for the the
species `spc`.
"""
function init_chrom(spc)
  # Find all chromophores in chemical formula
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

  return chrom
end #function init_chrom


"""
    init_ChromClass()

Initilise an empty array of length 19 with String arrays for each chromophore
class (to be appended by the MCM names of the species in each class) and return it.
"""
function init_ChromClass()

  # Define empty arrays for each chromophore class
  Ald = String[]; Ket = String[]; DiCar = String[]; Kete = String[];
  Nitro = String[]; DiNitro = String[]; Nit = String[]; DiNit = String[];
  PNit = String[]; DiPNit = String[]; PAN = String[];
  ROOH = String[]; DiROOH = String[]; PAA = String[]; SCI = String[];
  Poly = String[]; NoChr = String[]; Rad = String[]; New = String[]

  # Compile single arrays in overall array and return it
  return [Ald, Ket, DiCar, Kete, Nitro, DiNitro, Nit, DiNit, PNit, DiPNit, PAN,
          ROOH, DiROOH, PAA, SCI, Poly, NoChr, Rad, New]
end


"""
    group_CC(spc, mspc, chrom, CC)

From current species' GECKO-A name `spc` and the array with the number of chromophores
in the current molecule `chrome`, append the array with lists of species in each
chromophore class `CC` by the MCM name of the species `mspc`.

The function only deals with the current scenario and does not loop over a range
of scenarios.
"""
function group_CC(spc, mspc, chrom, CC)
  # Save MCM names in the respective chromophore classes
  # Exclude radicals
  if     spc[1:3] == "new"                                    push!(CC[19],mspc)
  # Exclude radicals, but not stabilised Criegee intermediates
  elseif length(matchall(r"\.",spc))>0 &&
         length(matchall(r"\.\(OO\.\)(?!\*)",spc))==0         push!(CC[18],mspc)
  elseif all(chrom.==0)                                       push!(CC[17],mspc)
  # Species with no chromophores
  elseif all(chrom[2:end].==0) && chrom[1]==1                 push!(CC[1], mspc)
  # Carbonyl and ketene classes
  elseif chrom[1]==0 && all(chrom[3:end].==0) && chrom[2]==1  push!(CC[2], mspc)
  elseif (chrom[1]>0 || chrom[2]>0) && all(chrom[3:end].==0)  push!(CC[3], mspc)
  elseif all(chrom[1:2].==0) && all(chrom[4:end].==0)         push!(CC[4], mspc)
  # Classes of nitrogen containing species
  elseif all(chrom[1:3].==0) && all(chrom[5:end].==0)
    if chrom[4]==1                                            push!(CC[5], mspc)
    else                                                      push!(CC[6], mspc)
    end
  elseif all(chrom[1:4].==0) && all(chrom[6:end].==0)
    if chrom[5]==1                                            push!(CC[7], mspc)
    else                                                      push!(CC[8], mspc)
    end
  elseif all(chrom[1:5].==0) && all(chrom[7:end].==0)
    if chrom[6]==1                                            push!(CC[9], mspc)
    else                                                      push!(CC[10],mspc)
    end
  elseif all(chrom[1:6].==0) && all(chrom[8:end].==0)         push!(CC[11],mspc)
  # Classes of species with a hydroperoxide group
  elseif all(chrom[1:7].==0) && all(chrom[9:end].==0)
    if chrom[8]==1                                            push!(CC[12],mspc)
    else                                                      push!(CC[13],mspc)
    end
  elseif all(chrom[1:8].==0) && all(chrom[9:end].==0)         push!(CC[14],mspc)
  # Stabilised Criegee intermediates
  elseif all(chrom[1:9].==0)                                  push!(CC[15],mspc)
  # Species with polyfunctional chromophores
  else                                                        push!(CC[16],mspc)
  end

  # Return the appended array with lists of species in each chromophore class
  return CC
end #function group_CC


"""
    init_OCclass()

Initilise an empty array of length 4 with String arrays for each O:C ratio class
(to be appended by the MCM names of the species in each class) and return it.
"""
function init_OCclass()
  oc1 = String[]; oc2 = String[]; oc3 = String[]; oc4 = String[]
  return [oc1,oc2,oc3,oc4]
end


"""
    group_OC(OCclasses,OCratio,gspc,mspc)

From current species' GECKO-A name `gspc` and the O:C ratio of the current
molecule `OCratio`, append the array with lists of species in each
O:C ratio class `OCclasses` by the MCM name of the species `mspc`.

The function only deals with the current scenario and does not loop over a range
of scenarios.
"""
function group_OC(OCclasses,OCratio,gspc,mspc)
  if searchindex(gspc,".") > 0
    # Exclude radicals
    return OCclasses
  elseif   gspc[1:3] == "new"
    # Add new generic species
    if     gspc[end] == '1'     push!(OCclasses[1],mspc)
    elseif gspc[end] == '2' push!(OCclasses[2],mspc)
    elseif gspc[end] == '3' push!(OCclasses[3],mspc)
    elseif gspc[end] == '4' push!(OCclasses[4],mspc)
    end
  # Assignment to O:C ratio classes for general cases
  elseif OCratio≤0.5  push!(OCclasses[1],mspc)
  elseif OCratio≤1.0  push!(OCclasses[2],mspc)
  elseif OCratio≤2.0  push!(OCclasses[3],mspc)
  else                push!(OCclasses[4],mspc)
  end

  # Return array with appended lists of O:C ratio classes
  return OCclasses
end



"""
    init_SizeClass()

Initilise an empty array of length 10 with String arrays for each size class
(to be appended by the MCM names of the species in each class) and return it.
"""
function init_SizeClass()
  sc1 = String[]; sc2 = String[]; sc3 = String[]; sc4 = String[]
  sc5 = String[]; sc6 = String[]; sc7 = String[]; sc8 = String[]
  sc9 = String[]; sc0 = String[]
  return [sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc0]
end


"""
    group_CN(SizeClasses,CN,gspc,mspc)

From current species' GECKO-A name `gspc` and the size (carbon number + ether
groups) of the current molecule `CN`, append the array with lists of species in each
size class `SizeClasses` by the MCM name of the species `mspc`.

The function only deals with the current scenario and does not loop over a range
of scenarios.
"""
function group_CN(SizeClasses,CN,gspc,mspc)
  if searchindex(gspc,".") > 0 || gspc[1:3] == "new"
    # Exclude radicals, Criegee intermediates, and new generic model species
    return SizeClasses
  # Assign size class according to number of carbon and ether groups
  elseif CN==1  push!(SizeClasses[1],mspc)
  elseif CN==2  push!(SizeClasses[2],mspc)
  elseif CN==3  push!(SizeClasses[3],mspc)
  elseif CN==4  push!(SizeClasses[4],mspc)
  elseif CN==5  push!(SizeClasses[5],mspc)
  elseif CN==6  push!(SizeClasses[6],mspc)
  elseif CN==7  push!(SizeClasses[7],mspc)
  elseif CN==8  push!(SizeClasses[8],mspc)
  elseif CN==9  push!(SizeClasses[9],mspc)
  elseif CN≥10  push!(SizeClasses[10],mspc)
  else
    # Warning, when failed
    println("Error determining size of $gspc ($mspc). Species ignored in size classes.")
  end

  # Return array with appended lists of size classes
  return SizeClasses
end


"""
    group_conc(output,class_data,classes)

Append dataframe `output` with species concentrations by concentrations of classes
with certain properties using the header names define in the array `classes`.
"""
function group_conc(output,class_data,classes)
  # Loop over chromophore classes
  for (i, class) in enumerate(classes)
    # Initilise concentration arrays for each chromophore class
    output[Symbol(class)] = zeros(length(output[1]))
    # Loop over species in each chromophore class
    for spc in class_data[i]
      # Add concentrations of individual species to chromophore class
      output[Symbol(class)] += output[Symbol(spc)]
    end
  end

  return output
end

end #module groupSPC
