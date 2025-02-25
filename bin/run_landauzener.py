import os
from shutil import copy2 as copy
#
import numpy as np
#
from pysurf.logger import get_logger
from pysurf.sampling import Sampling
from pysurf.setup import SetupBase
from pysurf.spp import SurfacePointProvider
from pysurf.utils import exists_and_isfile
#
from pysurf.utils.constants import fs2au
from pysurf.dynamics.base_propagator import PropagatorFactory
#
from colt import Colt, from_commandline


class Prop(Colt):

    _user_input = """
    # Total propagation time in fs
    time_final [fs] = 100 :: float

    # Time step in fs for the propagation
    timestep [fs] = 0.5 :: float

    # File with initial condition
    initial condition = init.db :: str

    # Number of total states
    n_states = :: int
    nghost_states = 0 :: int
    method = LandauZener :: str

    #Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: str

    restart = True :: bool

    additional = :: str, optional :: sts_mom, fosc
    """

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases("method", {name: method.colt_user_input
                                 for name, method in PropagatorFactory._methods.items()})

    def __init__(self, config):
        self.logger = get_logger('prop.log', 'prop')
        self.logger.header('PROPAGATION', config)
        sampling = Sampling.from_db('sampling.inp', config['initial condition'])

        self.nsteps = int(np.ceil(config['time_final [fs]']/config['timestep [fs]']))

        self.logger.info('Start propagation of trajectory')

        if config['additional'] is None:
            # get propagator
            propagator = PropagatorFactory._methods[config['method'].value](config['spp'],
                                                                            sampling,
                                                                            config['n_states'],
                                                                            nghost_states=config['nghost_states'],
                                                                            # properties = config['properties'],
                                                                            restart=config['restart'],
                                                                            logger=self.logger)

        else:
            propagator = PropagatorFactory._methods[config['method'].value](config['spp'],
                                                                            sampling,
                                                                            config['n_states'],
                                                                            nghost_states=config['nghost_states'],
                                                                            properties=[config['additional']],
                                                                            restart=config['restart'],
                                                                            logger=self.logger)
        propagator.run(self.nsteps, config['timestep [fs]']*fs2au)

    @classmethod
    def from_config(cls, config):
        return cls(config)


class SetupPropagation(SetupBase):

    folder = 'prop'
    subfolder = 'traj'

    _user_input = """
    # Number of trajectories for the propagation
    n_traj = -1 :: int

    # sampling input
    sampling_input = sampling.inp :: existing_file

    # Database containing all the initial conditions
    sampling_db = sampling.db :: existing_file

    # Filepath for the inputfile of the Surface Point Provider
    spp = spp.inp :: file

    # initial excited state for the trajectory
    initial state =  :: int

    # Filepath for the inputfile of the Propagation
    prop = prop.inp :: file

    # Decide whether database for the propagation should be copied to the trajectory folder
    copy_db = none :: str

    ## Run LZ or FSSH
    #tsh_method =  :: str :: FSSH, LZ
    """

    def __init__(self, config):
        #self.tsh_method = config['tsh_method']
        """ Class to create initial conditions due to user input. Initial conditions are saved
            in a file for further usage.
        """
        logger = get_logger('setup_landauzener.log', 'setup_landauzener')
        SetupBase.__init__(self, logger)


        # Open DB of initial conditions once, so that it is available
        sampling = Sampling.from_db(config['sampling_input'], config['sampling_db'], logger=logger)

        self.moldenfile = None
        sampling_config = Sampling.generate_user_input(config['sampling_input']).check_only()
        if sampling_config['method'] == 'Wigner':
            if sampling_config['method']['from'] == 'molden':
                self.moldenfile = sampling_config['method']['from']['moldenfile']

        # Make sure that inputfile for the SPP exists and is complete

        if exists_and_isfile(config['spp']):
            lconfig = config['spp']
        else:
            lconfig = None
        spp_config = SurfacePointProvider.generate_input(config['spp'], config=lconfig)
        if spp_config['mode'] == 'ab-initio':
            self.model = None
            self.basis = spp_config['mode']['software']['basis']
        elif spp_config['mode']['model'] == 'LVC':
            self.model = "LVC"
            self.basis = None
        else:
            self.basis = None
            self.model = None


        spp_config = SurfacePointProvider.generate_input(config['spp'], config=lconfig)

        # Make sure that inputfile for RunTrajectory exists and is complete
        if exists_and_isfile(config['prop']):
            lconfig = config['prop']
        else:
            lconfig = None

        prop = Prop.generate_input(config['prop'], config=lconfig)
        if config['n_traj'] == -1:
            ntraj = len(sampling._db)
        else:
            ntraj = config['n_traj']
        if sampling.nconditions < ntraj:
            logger.error(f"Too few initial conditions in {config['sampling_db']}")

        self.setup_folders(range(ntraj), config, sampling, prop['initial condition'])



    @classmethod
    def from_config(cls, config):
        return cls(config)

    def setup_folder(self, number, foldername, config, sampling, initdb):
        copy(config['prop'], foldername)
        copy(config['spp'], foldername)
        if self.model == 'LVC':
            copy('pyrmod6.inp', foldername)

        if self.basis == 'gen':
            copy('basis_gen.ini', foldername)

        if config['copy_db'] != 'none':
            copy(config['copy_db'], foldername)

        copy(config['sampling_input'], foldername)
        if self.moldenfile is not None:
            copy(self.moldenfile, foldername)
        #setup new database
        initname = os.path.join(foldername, initdb)
        new_sampling = Sampling.create_db(config['sampling_input'], initname, sampling.info['variables'], sampling.info['dimensions'], sampling.system, sampling.modes, model=sampling.model, sp=False)
        #copy condition to new db
        condition = sampling.get_condition(number)
        new_sampling.write_condition(condition, 0)
        new_sampling.set('currstate', config['initial state'], 0)


@from_commandline("""
# Which option to select
variant = :: str 
[variant(run)]
[variant(setup)]
""")
def main(variant):
    if variant == 'run':
        run = Prop.from_questions(config='prop.inp')
    elif variant == 'setup':
        SetupPropagation.from_questions(config="setup_landauzener.inp")
    else:
        raise ValueError("do not understand variant")


if __name__ == '__main__':
    main()
