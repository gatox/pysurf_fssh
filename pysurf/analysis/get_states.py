#! /data/ehrmaier/anaconda3/bin/python3
import sys
import os 
import numpy as np
import click

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.utils.chemutils import get_atom_from_mass

from pysurf.utils.constants import au2ev



def write_state(state, step):
    string = ''
    string += 'step {0} : {1}\n'.format(step, state)
    return string


@click.command()
@click.option('-f', 'infile', default='prop.db')
@click.option('-o', 'outfile', default='states.dat')
def get_states_command(infile, outfile):
    get_states(infile, outfile)

def get_states(infile, outfile):
    if not(os.path.isfile(infile)):
        print('Error: infile path does not exist! ' + infile)
        exit()
    
    db = Database.load_db(infile)
    
    with open(outfile, 'w') as output:
        step = 0
        for state in db['curr_state']:
            output.write(write_state(state[0], step))
            step += 1
    
if __name__=="__main__":
    get_states_command()