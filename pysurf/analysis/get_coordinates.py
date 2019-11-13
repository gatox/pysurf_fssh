#! /data/ehrmaier/anaconda3/bin/python3
import sys
import os 
import cli
import numpy as np
import click

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.utils.chemutils import get_atom_from_mass

from pysurf.utils.constants import bohr2angstrom

def au2ang(x):
    return x * bohr2angstrom

def write_xyz(atoms, coord, step):
    string = str(len(coord)) + '\n'
    string += 'step {0} \n'.format(step)
    for i in range(len(coord)):
        string += '{0:s}  {1:12.8f}  {2:12.8f}  {3:12.8f}\n'.format(atoms[i], *au2ang(coord[i]))
    return string

def write_coord(coord, step):
    string = str(len(coord)) + '\n'
    string += 'step {0} \n'.format(step)
    np.vectorize(str)
    string += np.array2string(coord, separator=',   ', precision=5).strip(']').strip('[')
    string += '\n'
    string += '\n'
    return string


@click.command()
@click.option('-f', 'infile', default='prop.db')
@click.option('-o', 'outfile', default='coord.xyz')
def get_coordinates_command(infile, outfile):
    get_coordinates(infile, outfile)

def get_coordinates(infile, outfile):
    if not(os.path.isfile(infile)):
        print('Error: infile path does not exist! ' + infile)
        exit()
    
    db = Database.load_db(infile)
    try:
        mass = db['mass']
        if len(mass.shape) == 1:
            model = True
        else:
            model = False
    except:
        print('Masses could not be found in DB!')
        mass = None
        model = True
    
    atoms=[]
    if model is True:
        for m in range(len(db['coord'][0])):
            atoms+=['Q']
    if model is False:
        for m in mass[:,0]:
            atoms += [get_atom_from_mass(m)]
    
    
    with open(outfile, 'w') as output:
        step = 0
        for coord in db['coord']:
            if model is False:
                output.write(write_xyz(atoms, coord, step))
            else:
                output.write(write_coord(coord, step))
            step += 1
    
if __name__=='__main__':
    get_coordinates_command()