from pysurf.logger import get_logger 
from pysurf.sampling import Sampling

from pysurf.utils import exists_and_isfile
from pysurf.spp import SurfacePointProvider
from colt import Colt, from_commandline


class SinglePointCalculation(Colt):
    """ The SinglePointCalculation class asks the SPP for information at one specific geometry.
        
        It takes its information from the init_db and writes the results to it.
        It was designed specifically for single point calculations for the spectrum
        Properties:
        -------
        Classmethods:
            
            from_inputfile:
                classmethod to initialize the class from an inputfile
            from_config:
                classmethod to initialize the class, e.g. used in cls.from_commandline
            
        -------
        Methods:
    """

    _user_input = """
    # Number of states
    nstates =  :: int
    nghost_states = 0 :: int
    # Database containing the geometry
    init_db = init.db :: str
    # Properties that should be calculated
    properties = [energy, fosc] :: list
    # SPP inputfile
    spp = spp.inp :: file
    # sampling inp
    sampling_input = sampling.inp :: existing_file
    """

    def __init__(self, config, logger=None):
        """ 
            Args:
                config, ColtObj:
                    The config contains the information from the colt 
                    questions of the class
                logger, Logger:
                    Logger for the class. If not provided a new logger is created.
        """

        if logger is None:
            self.logger = get_logger('sp_calc.log', 'sp_calc')
            self.logger.header('Single Point Calculation', config)
        else:
            self.logger = logger
        #
        self.logger.info(f"Taking information from {config['init_db']}")
        sampling = Sampling.from_db(config['sampling_input'], config['init_db'], logger=self.logger, sp=True)

        if not(exists_and_isfile(config['spp'])):
            spp_config = SurfacePointProvider.generate_input(config['spp'])
        else:
            spp_config = SurfacePointProvider.generate_input(config['spp'], config=config['spp'])
            
   
        self.logger.debug(f"Setting up SPP with {config['spp']}")
        if sampling.molecule is not None:
            spp = SurfacePointProvider.from_config(spp_config, 
                                      config['properties'],
                                      config['nstates'],
                                      sampling.natoms,
                                      nghost_states=config['nghost_states'],
                                      atomids = sampling.atomids)
        else:#model calculation
            spp = SurfacePointProvider,from_config(spp_config, 
                                      config['properties'],
                                      config['nstates'],
                                      sampling.nmodes)

        crd = sampling.get_condition(0).crd

        #check that DB can take the results
        self.logger.debug(f"Checking whether {config['init_db']} can take the results")
        self._check_res_db(config, sampling)
        
        with self.logger.info_block("SPP calculation"):
            res = spp.request(crd, config['properties'], states=[st for st in range(config['nstates'])])

        self.logger.info(f"Writing results to: {config['init_db']}")
        for prop, value in res.iter_data():
            sampling.set(prop, value)
            
    @classmethod
    def from_config(cls, config):
        return cls(config)

    @classmethod
    def from_inputfile(cls, inputfile):
        config = cls.generate_input(inputfile, config=inputfile)
#        config['inputfile'] = inputfile
        logger = get_logger('sp_calc.log', 'sp_calc')
        logger.header('Single Point Calculation', config)
        return cls(config, logger)

    def _check_res_db(self, config, sampling):
        if exists_and_isfile(config['init_db']):
            info = sampling.info
            check1 = all(item in info['variables'] for item in config['properties'])
            if 'nstates' in info['dimensions']:
                if (info['dimensions']['nstates'] >= config['nstates']): 
                    if 'natoms' in info['dimensions']:
                        if(info['dimensions']['natoms'] == sampling._db.info['dimensions']['natoms']):
                            check2 = True
                        else:
                            check2 = False
                    if 'nmodes' in info['dimensions']:
                        if(info['dimensions']['nmodes'] == sampling._db.info['dimensions']['nmodes']):
                            check2 = True
                        else:
                            check2 = False
                else:
                    check2 = False
            else:
                check2 = False
            if check1 and check2:
                return 
            else:
                self.logger.error(f"Given database is not appropriate for the results: {config['init_db']}")
                    
@from_commandline("""
inputfile = sp_calc.inp :: file 
""")
def command_sp_calc(inputfile):
    """ Setting up initial conditions according to the inputfile.
    If inputfile doesn't exist, colt will ask all necessary questions
    """
    sp_calc = SinglePointCalculation.from_inputfile(inputfile)
#    sampling = Sampling.from_db('sampling.db')


if __name__=="__main__":
    command_sp_calc()
