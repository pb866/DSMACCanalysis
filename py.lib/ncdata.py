import netCDF4,re
from netCDF4 import Dataset

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
