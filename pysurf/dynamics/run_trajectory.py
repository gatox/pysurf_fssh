import os
import numpy as np

from collections import namedtuple

from ..logger import Logger, get_logger
from ..database.database import Database
from ..database.dbtools import DBVariable
from ..utils.constants import fs2au
from ..utils.strutils import split_str
from ..spp import SurfacePointProvider
from .base_propagator import PropagatorFactory
from ..sampling import Sampling

#
from colt import Colt


class RunTrajectory(Colt):
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

#    properties = energy, gradient :: list
    """

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases(
            "method",
            {
                name: method.colt_user_input
                for name, method in PropagatorFactory._methods.items()
            },
        )

    def __init__(self, config):
        self.logger = get_logger("prop.log", "prop")
        self.logger.header("PROPAGATION", config)
        sampling = Sampling.from_db(config["initial condition"])

        self.nsteps = int(np.ceil(config["time_final [fs]"] / config["timestep [fs]"]))

        self.logger.info("Start propagation of trajectory")

        # get propagator
        propagator = PropagatorFactory._methods[config["method"].value](
            config["spp"],
            sampling,
            config["n_states"],
            nghost_states=config["nghost_states"],
            # properties = config['properties'],
            restart=config["restart"],
            logger=self.logger,
        )
        propagator.run(self.nsteps, config["timestep [fs]"] * fs2au)

    @classmethod
    def from_inputfile(cls, inputfile):
        quests = cls.generate_user_input(config=inputfile)
        config = quests.ask(inputfile)
        return cls(config)
