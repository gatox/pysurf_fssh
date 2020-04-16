import os
import numpy as np
from shutil import copy2

from pysurf.colt import Colt
from pysurf.colt import FromCommandline
from pysurf.utils.osutils import exists_and_isfile
from pysurf.logger import get_logger
from pysurf.sampling import Sampling

class SetupPropagation(Colt):
    _questions = """
    # Total propagation time in fs
    propagation time in fs = 100 :: float

    # Time step in fs for the propagation
    time step in fs = 0.5 :: float
    
    # Number of trajectories for the propagation
    n_traj = 100 :: int

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file

    # Decide whether database for the propagation should be copied to the trajectory folder
    copy_db = yes :: str :: [yes, no]

    [copy_db(yes)]
    db_file = none :: str
    """
    

    def __init__(self, inputfile):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        self.inputfile = inputfile
        self.propagationfolder = 'prop'
        self.trajfolder = 'traj_'
        self.logger = get_logger('setup_propagation.log', 'setup_propagation')
        
        quests = self.generate_questions("SURFACE HOPPING", config=inputfile)
        self.config = quests.ask(inputfile)

        if not(os.path.isdir(self.propagationfolder)):
            os.mkdir(self.propagationfolder)
        
        # Open DB of initial conditions once, so that it is available
        self.sampling = Sampling.from_db(self.config['sampling_db'])

        for i in range(self.config['n_traj']):
            self._setup_trajectory_folder(i)

    def _setup_trajectory_folder(self, number):
        foldername = os.path.join(self.propagationfolder, self.trajfolder + '{:04d}'.format(number))
        if os.path.isdir(foldername):
            self.logger.info('Folder {0} already exists. Skipping trajectory'.format(foldername))
            return

        os.mkdir(foldername)

        copy2(self.inputfile, foldername)
        copy2(self.config['spp'], foldername)

        if self.config['copy_db'] == 'yes':
            copy2(self.config['copy_db']['db_file'], foldername)

        initname = os.path.join(foldername, 'init_{:04d}.db'.format(number))
        self.sampling.export_condition(initname, number)


@FromCommandline("""
inputfile = propagation.inp :: file
""")
def command_setup_propagation(inputfile):
    SetupPropagation(inputfile)

if __name__=="__main__":
    command_setup_propagation()
