import netCDF4,re
from netCDF4 import Dataset
import pandas as pd

def get(filename):
    nc = Dataset(filename,'r')
    print nc.date, '\n', nc.description,'\ngroup 0 selected'
    # print 'Select Simulation: \n\n'
    for i,g in enumerate(nc.groups): print i , ' - ', g
    group = tuple(nc.groups)[0]#[int(input('Enter Number \n'))]
    print group, 'took', nc.groups[group].WALL_time, 'seconds to compute.'

    specs = nc.groups[group].variables['Spec'][:]
    specs_columns = nc.groups[group].variables['Spec'].head.split(',')
    rates = nc.groups[group].variables['Rate'][:]
    rates_columns = nc.groups[group].variables['Rate'].head.split(',')

    nc.close()

    # di= dict([[specs_columns[i],i] for i in xrange(len(specs_columns))])

    return {'specs':specs,'rates':rates,
    'sc':specs_columns,'rc':rates_columns}


def load(filename):
    # Load data from file and print description
    nc = Dataset(filename,'r')
    group = tuple(nc.groups)[0]
    print group, 'took', nc.groups[group].WALL_time, 'seconds to compute.'
    print nc.date, '\n', nc.description

    # Save concentrations and rates to Pandas DataFrames
    specs = pd.DataFrame(nc.groups[group].variables['Spec'][:])
    specs.columns = nc.groups[group].variables['Spec'].head.split(',')
    rates = pd.DataFrame(nc.groups[group].variables['Rate'][:])
    rates.columns = nc.groups[group].variables['Rate'].head.split(',')[:]

    # Correct first time step and use time as index
    specs.TIME[0] = 2*specs.TIME[1] - specs.TIME[2]
    rates.TIME[0] = 2*rates.TIME[1] - rates.TIME[2]
    specs.index = pd.to_datetime(specs.TIME, unit='s')
    rates.index = pd.to_datetime(specs.TIME, unit='s')

    # Use DateTime format for time
    t = specs.TIME.index
    specs.TIME = t
    rates.TIME = t

    return specs, rates
