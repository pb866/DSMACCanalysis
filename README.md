Module DSMACCanalysis
=====================

Julia script to analyse and display DSMACC model output.

Content
-------

- Module __*ropa\_tool*__ to analyse chemical source and sink fluxes of DSMACC output
- Includes:
  - Module _FluxAnalysis_ to determine chemical source and sink fluxes for each
    species
  - Module *prepare\_ropa* to arrange DSMACC output from different scenarios in
    an array of data frames and analyse each reaction finding the educts and
    products with the corresponding stoichiometric indices
  - Module _NCload_ to load species concentrations and reaction rates into a
    dictionary with data frames
  - Function *test\_file* from module _fhandle_ in repository
    [auxdata](https://github.com/pb866/auxdata.git) to check existance of an nc file
  - Julia modules _DataFrames_ and _PyCall_
  - Function _input_ from Julia module _Juno_ to read user input
- Module __*groupSPC*__ to lump the NMVOC organic pool by certain properties
- Properties currently available:
  - Molecule's size (carbon and ether group skeleton)
  - O:C ratio
  - Chromophore class
- Includes:
  - Function _rdinp_ from module _fhandle_ in repository
    [auxdata](https://github.com/pb866/auxdata.git) to read in data files
  - Julia modules _DataFrames_
  - Function _input_ from Julia module _Juno_
- Module __*DSMACCplot*__ to plot DSMACC model results from specifications of an
  input file
- Includes:
  - Module _groupSPC_ to lump species according to their O:C ratio, their
    molecule's size, and their chromophore class
  - Module _ropa\_tool_ to analyse source and sink fluxes of each species
  - Module _jlplot_ to compile plots of DSMACC model output in a pdf
  - Module _plot\_ropa_ to plot time-resolved chemical sink and source fluxes
    for a species from the ROPA analysis
  - Module _make\_plots_ from repository
    [auxdata](https://github.com/pb866/auxdata.git) to generate python plots
    (with PyPlot)
  - Module _NCload_ to load species concentrations and reaction rates into a
    dictionary with data frames
  - Module _fhandle_ in repository [auxdata](https://github.com/pb866/auxdata.git)
    for input file handling
  - Julia modules DataFrames, PyCall, and PyPlot
  - Function _input_ from Julia module _Juno_


Installation
------------

The script is written for and test with
[Julia version 0.6.2](https://julialang.org/downloads/). If not installed, install Julia. Using `Pkg.add("<package name>")` the following packages need to be installed
in Julia:

- PyCall
- PyPlot
- DataFrames
- Juno

The script will analyse DSMACC model output from the
[DSMACC-testing version](https://github.com/pb866/DSMACC-testing.git).
Create a git submodule in the `AnalysisTools` folder named `DSMACCanalysis`
(or clone this repository into the AnalysisTools folder).

You will furthermore need the modules `fhandle` and `make_plots` from the
[auxdata](https://github.com/pb866/auxdata.git) repository. Clone the repository
to a directory of your liking and specify the directory paths in the preambles of
the following modules:

- DSMACCplot (ll. 27 – 38)
- jl.mod/groupSPC (ll. 46 – 60)
- jl.mod/jlplot.jl (ll. 40 – 52)
- jl.mod/NCload (ll. 31 – 43)
- jl.mod/ropa/plot_ropa (ll. 33 – 44)


Either adjust the folder path in one of the if statements to the directory of the
repository auxdata (3 times) as

```
directory/to/repo/jl.mod
```

or add a new case.


Usage of the script
-------------------

### The Julia ROPA tool

Function `ropa` in module `ropa_tool` (in the jl.mod/ropa folder) is used to
analyse the sink and source fluxes for every species in a series of model runs of
specified DSMACC model output netCDF files. It can be used as a standalone tool
to return arrays with DataFrames for each scenario holding the source and sink
fluxes for each species in the specified mechanisms and the species concentrations.
Call the ropa tool with:

```julia
source, sink, specs = function ropa(scenarios=[]; specs=[], rates=[], cycles="reduce")
```

Where the argument scenarios is a string list with the species names (and folder
positions). If you already retrieved the species concentrations and reaction rates
as array of DataFrames for each scenario, you can omit the scenarios argument and
specify the keyword arguments `specs` for the species concentrations and `rates`
for the reaction rates instead. Furthermore, the ROPA tool calculates net fluxes
of inorganic NOx and Ox recycling as those fluxes would otherwise overlay the true
sinks and sources for important Ox and NOx species. If you don't want net fluxes,
set the keyword argument `cycles` to `"full"`. Otherwise, the following reactions
will be shown as net main sources and sinks for _O<sub>3</sub>, O(<sup>1</sup>D),
O(<sup>3</sup>P), NO, and NO<sub>2</sub>_:

- O<sub>3</sub> ⟶ O(<sup>1</sup>D)
- O(<sup>1</sup>D) ⟶ O(<sup>3</sup>P)
- O<sub>3</sub> ⟶ O(<sup>3</sup>P)
- O(<sup>3</sup>P) ⟶ O<sub>3</sub>
- NO<sub>2</sub> ⟶ NO + O(<sup>3</sup>P)
- NO + O<sub>3</sub> ⟶ NO<sub>2</sub>

The output data is organised as arrays of type `Any` holding data frames of each
scenario. The `concs` DataFrames hold columns for the concentrations at every time
step of every species using the DSMACC species name as column header.

The `sinks` and `sources` DataFrames hold an n×2 matrix with the all the species
names in the first column. The second column holds another n×2 matrix with the
DSMACC reaction (from the labels of the `rates` DataFrames) in the first column
and the chemical flux data in the second column. For each species all respective
sink and source reaction fluxes are listed in the matrices in the second column
of the outer matrix.


### The Julia module _groupSPC_

The module is used by the plotting script _DSMACCplot.jl_ (see below). It groups
the pool of NMVOC species according to certain properties. The lumped concentrations
are added to the DataFrames with the species concentrations in each scenario using
the column headers given below. These column headers can be used in the input file
with the plot specifications to plot the concentrations of the respective property.

The following properties exist:


### Molecule's size

Molecules are lumped by size. The size is determined by the carbon number plus the
number of ether groups in a molecule (making up the molecule's skeleton). 10 size
classes exist, with class 1 to 9 holding that number of carbon + ether groups and
class 10 holding 10 or more carbon + ether groups. The keyword for each class is
`sc#`, where `#` is `1` to `9` for the first 9 classes and `0` for class 10.


### O:C ratio

Molecules are grouped by O:C ratio into 4 classes with increasing O:C ratios and,
hence, expected increased solubility. In the first class, pure gas phase species
are expected, the last class exists of species, which are expected dominantly or
exclusively in the aqueous phase. Classes 2 and 3 exist of semi-volatile species
with increasing solubility.

- `oc1`:       O:C ≤ 0.5
- `oc2`: 0.5 < O:C ≤ 1.0
- `oc3`: 1.0 < O:C ≤ 2.0
- `oc4`: 2.0 < O:C


### Chromophore classes

Stable (non-radical) species are furthermore lumped by the chromophores they hold.
Classes exist for species holding a single chromophore or chromophores of only 1
type. Another class exists for species with several chromophores of different type
(polyfunctional chromophores) or no chromophores at all. Further special cases are
inorganic species, organic radical species, and generic 'new' model species from
test scenarios with a new photolysis protocol of the Master Chemical Mechanism
([MCM](http://mcm.leeds.ac.uk/MCMv3.3.1/home.htt)).

The following classes exist:

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


### The Julia plotting script _DSMACCplot.jl_

Although the previous modules can be used as standalone versions, they are
conveniently incorporated in the _DSMACCplot.jl_ script, where their output is
displayed graphically. Standalone versions of the above modules are only needed,
if further data manipulation is desired. To display DSMACC model output, the
following plot types exist in the current version:

- Line plots for species concentrations and reaction rates
- Stacked area plots for species concentrations (meant for lumped sum species)
  for up to 2 different scenarios using border lines with different line styles
  and different colour schemes for different scenarios
- Stacked area plots for time-resolved sink and source fluxes for a particular
  species in a chosen scenario

Run the script using

```
julia DSMACCplot.jl [[<path/>]<input file>] [<path/>[<output scenario name>]]
```

where the first script argument defines the name (and directory) of a text file
with specifications of the desired plots. The second arguments defines the file
name (without the file ending) and folder path of the output pdf with the desired
plots. If the input file name is missing in the first argument or the file doesn't
exist, you will be asked for user input during the execution of the script.

A pdf with all the specified plots will be produced in the directory from which the
script was called (i.e., either the DSMACC main folder or
./AnalysisTools/DSMACCanalysis) unless a different folder path is specified in
the second script argument. If no name was specified in the second script argument,
the labels of each scenario joined by an underscore (`_`) are used as file name.


The input file with plot specifications
---------------------------------------

The input file is devided in 4 sections, each addressed by a section caption. It
is important to keep the order of the sections the following:

1. `Scenarios` section
2. `Settings` section (optional)
3. `Plotting` section
4. `Comments` section (optional)


### Scenarios section

The first section is labelled with the keyword `Scenarios:`. On the first line after
the caption, list all the netCDF files, you want to consider for plotting. If the
files are saved in the folder `./save/results` you can omit the folder path, otherwise
folder paths are mandatory. List all files on one line. The following separators are
allowed on horizontal lists on the same line throughout the input file:

- whitespace (spaces, tabs)
- commas (`,`)
- semicolons (`;`)

On the second line you can specify labels for each scenario, which are used in the
figure legends. If this line is obsolete, file names are used as labels (without the
file ending). Use those automatic labels to address, which scenarios you want to plot
(see [plotting section](#plotting-section)).


### Settings section

The second section is optional and is introduced by the keyword `Settings`. It lists
general specifications of the plots. Currently a lower and upper cut-off can be
defined. Chemical fluxes smaller than the lower cut-off are combined in the plots.
If fluxes larger than the upper cut-off exist, a second zoom plot will be created,
where fluxes above the upper cut-off are omitted to allow the view of the smaller
fluxes. Define the parameters using the keyword `cut-off:` followed by a list of
the lower and upper threshold on that line. If you don't want this feature, set the
lower cut-off to `0.0` and the upper cut-off to `1.0`. If you don't specify thresholds,
standard values of `0.05` and `0.7` are used for the lower and upper cut-off,
respectively.

Furthermore, net fluxes of inorganic main NOx and Ox cycles can be formed during
the ROPA analysis as described above. By default, this feature is set. To illustrate
this in the input file use the keyword `cycles:` in the settings section followed
by `reduce`. You may omit this line, if you want net fluxes. Use any other phrase
after the keyword `cycles:` such as `full` to shut off this feature.


### Plotting section

The actual plot specifications are given in the third section started by the caption
`Plotting:`. Different plot types and plots can be specified. Plot types are
introduced by a caption line, specifying for which scenarios plots are to be produced,
what type of plot, and the units used in the plots. On the following lines, species
or reactions can be specified, for which the data should be plotted. Different plot
types can be specified in one file and must be separated by at least one blank line,
while blank lines are prohibited within one plot type section.

The following plots are currently available and can be specified as follows:


#### Line plots

For species concentrations use the heading:

    <scen1> ... <scen n>: specs/unit


For reaction rates use the heading:

    <scen1> ... <scen n>: rates/unit

Use a horizontal list of scenarios, if you want to plot concentrations of the same
species for more than one scenario and the respective keyword for concentrations or
rates after a colon (`:`). Use a slash after the keyword and a keyword for the unit
to define units (see below).

Specify the species or rates to be plotted using the DSMACC labels from the netCDF
file for the species and rates. For every plot you want, start a new line and specify
all species that go into the plot in a vertical list (see also example file plot.inp).


#### Stacked area plots of concentrations

To represent the total NMVOC mass differentiated by certain properties, stacked
area plots exist. For this plot type, a transparent area plot with opaque boundary
lines is used. This plot type can be defined for up to two different scenarios.
Different line styles of the boundaries and different colour schemes will be used
for the different scenarios. Specify the plot type by the heading

    <scen1>[, <scen 2>]: stack/unit

Again, list data specification for different plots on different lines and the
(lumped) species for each plots as a horizontal list on 1 line.


#### Time-resolved ROPA plots

You can plot time-resolved sink and source fluxes from the ROPA analysis in stacked
area plots, where sink fluxes are negative and source fluxes are positive. Start
this plot type by the head line

    <scen1>, ... <scen n>: fluxes

followed by a list of species. The list of species can be on separate lines or in
a horizontal list on a single lines or mixtures of both.

In any case, flux plots are plotted in separate plots for every species and every
scenario no matter of the format in the input file with the plot specifications.


#### Unit specifications

Units can be specified after the plot keyword using a slash and the unit keyword.
The below unit keywords/units are currently allowed. The standard unit is
molecules cm<sup>-3</sup> for concentrations, cm<sup>3</sup> molecule<sup>-1</sup>
s<sup>-1</sup> for rates, and molecules cm<sup>-3</sup> s<sup>-1</sup> for fluxes.
For standard units the keyword can be omitted in the plot specifications. Flux plots
currently only exist in the default unit.

- `mlc` or `cm-3` (default): molecules cm<sup>-3</sup> or cm<sup>3</sup>
  molecule<sup>-1</sup> s<sup>-1</sup> or molecules cm<sup>-3</sup> s<sup>-1</sup>
- `ppm`: ppm or ppm<sup>1-n</sup> h<sup>-1</sup> (n … order of reaction)
- `ppb`: ppb or ppb<sup>1-n</sup> h<sup>-1</sup>
- `ppt`: ppt or ppt<sup>1-n</sup> h<sup>-1</sup>



### Comments section

The last section is for comments only and is ignored by the actual plotting script.
It is optional and can be started by the caption `Comments:`.



Version history
===============

Version 0.6.1
-------------

- Bugfix: Plotting previously omitted last line of species in the last plot type
  when Comments section is missing
- Bugfix: Correcting paths for module auxdata to be recognised on earth0

Version 0.6
-----------

- Julia ROPA tool for the time-resolved analysis of chemical sink and source fluxes
  for each species in a list of defined scenarios
- Julia module `groupSPC` to lump concentrations of NMVOC by the properties
  molecule's size, O:C ratio, and chromophore class
- Julia script `DSMACCplot.jl` to generate a pdf with line plots of species
  concentrations or reaction rates, stacked area plots (with boundary lines)
  of species concentrations (especially lumped concentrations from groupSPC),
  and time-resolved sink and source flux plots from the ROPA analysis
- Unit conversions from molecules cm<sup>-3</sup> / cm<sup>3</sup>
  molecule<sup>-1</sup> s<sup>-1</sup> to ppm / ppm h<sup>-1</sup>,
  ppb / ppb h<sup>-1</sup>, and ppt / ppt h<sup>-1</sup> are available for line
  plots and stacked area plots of concentrations
