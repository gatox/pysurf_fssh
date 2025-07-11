import numpy as np

from .sampling_db import SamplingDB
from .base_sampler import SamplerFactory

from ..utils import exists_and_isfile
from ..logger import Logger, get_logger
from ..system import Molecule

#
from colt import Colt


class Sampling(Colt):
    """Sampling is the header class for the sampling routines. It asks the main questions, selects
    the sampler and reads and writes the conditions to the sampling database.
    """

    _user_input = """
    # Maximium Number of Samples
    n_conditions_max = 100 :: int

    method = :: str

    # Database containing all the initial conditions
    sampling_db = sampling.db :: file

    """

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases(
            "method",
            {
                name: method.colt_user_input
                for name, method in SamplerFactory._methods.items()
            },
        )

    @classmethod
    def from_config(cls, config, logger=None, fillit=False):
        if exists_and_isfile(config["sampling_db"]):
            logger.info(f"Found existing database {config['sampling_db']}")
            db = SamplingDB.from_db(config["sampling_db"])
            logger.info(
                f"There are already {db.nconditions} conditions in the database"
            )
            sampler = None
        else:
            logger.info(f"Loading sampler {config['method'].value}")
            sampler = cls._get_sampler(config["method"], start=0)
            logger.info(f"Creating new database {config['sampling_db']}")
            db = SamplingDB.from_sampler(config, sampler)
        return cls(
            config, db, db.dynsampling, sampler=sampler, logger=logger, fillit=fillit
        )

    @classmethod
    def from_inputfile(cls, inputfile, fillit=True):
        logger = get_logger("sampling.log", "sampling")
        # Generate the config
        if exists_and_isfile(inputfile):
            config = cls.generate_input(inputfile, config=inputfile)
        else:
            config = cls.generate_input(inputfile)
        logger.header("SAMPLING", config)
        logger.info(f"Took information from inputfile {inputfile}")
        return cls.from_config(config, logger=logger, fillit=fillit)

    @classmethod
    def create_db(
        cls,
        sampling_config,
        dbfilename,
        variables,
        dimensions,
        system,
        modes,
        model=False,
        sp=False,
        logger=None,
        fillit=False,
    ):
        db = SamplingDB.create_db(
            dbfilename,
            variables,
            dimensions=dimensions,
            system=system,
            modes=modes,
            model=model,
            sp=sp,
        )
        # config = db.get_config()
        sampling_config = cls.generate_user_input(sampling_config).check_only()
        sampling_config.update({"sampling_db": dbfilename})
        return cls(
            sampling_config, db, db.dynsampling, logger=logger, sp=sp, fillit=fillit
        )

    @classmethod
    def from_db(cls, sampling_config, dbfilename, logger=None, sp=False, fillit=False):
        db = SamplingDB.from_db(dbfilename)
        # config = db.get_config()
        sampling_config = cls.generate_user_input(sampling_config).check_only()
        sampling_config.update({"sampling_db": dbfilename})
        return cls(
            sampling_config, db, db.dynsampling, logger=logger, sp=sp, fillit=fillit
        )

    def __init__(
        self, config, db, dynsampling, sampler=None, logger=None, sp=False, fillit=False
    ):
        """Sampling always goes with a database, if not needed use Sampler class"""
        self._db = db
        if logger is None:
            self.logger = get_logger("sampling.log", "sampling")
            self.logger.header("SAMPLING", config)
        else:
            self.logger = logger

        if sampler is None:
            self.logger.info(f"Loading sampler {config['method'].value}")
            self.sampler = self._get_sampler(config["method"], start=self.nconditions)
        else:
            self.sampler = sampler

        if "n_conditions_max" not in config:
            n_conditions_max = config["n_conditions"]
        else:
            n_conditions_max = config["n_conditions_max"]

        if sp is True:
            n_conditions = self.sampler.get_number_of_conditions(1)
            self.write_conditions(1)
        else:
            n_conditions = self.sampler.get_number_of_conditions(n_conditions_max)

            # check for conditions
            if fillit is True:
                if self.nconditions < n_conditions:
                    # setup sampler
                    self.logger.info(
                        f"Adding {n_conditions-self.nconditions} additional entries to the database"
                    )
                    self.add_conditions(n_conditions - self.nconditions)

    def add_conditions(self, nconditions, state=0):
        # TODO
        # Take the random seed and random counter from the database to
        # assure the consistency with all
        # the previous conditions
        for _ in range(nconditions):
            cond = self.sampler.get_condition()
            if cond is None:
                self.logger.error("Sampler has no more conditions")
            self._db.append_condition(cond)

    def write_conditions(self, nconditions, state=0):
        # TODO
        # Take the random seed and random counter from the database to
        # assure the consistency with all
        # the previous conditions
        for _ in range(nconditions):
            cond = self.sampler.get_condition()
            if cond is None:
                self.logger.error("Sampler has no more conditions")
            self._db.write_condition(cond, 0)

    def write_condition(self, condition, idx):
        self._db.write_condition(condition, idx)

    def append_condition(self, condition, idx):
        self._db.append_condition(condition, idx)

    def append(self, key, value):
        self._db.append(key, value)

    def get(self, key, idx=None):
        return self._db.get(key, idx)

    def set(self, key, value, idx=None):
        self._db.set(key, value, idx)

    def __iter__(self):
        self._start = 1  # skip equilibrium structure
        return self

    def __next__(self):
        cond = self.get_condition(self._start)
        if cond is not None:
            self._start += 1
            return cond
        raise StopIteration

    def get_condition(self, idx):
        return self._db.get_condition(idx)

    @staticmethod
    def _get_sampler(config, start=0):
        return SamplerFactory._methods[config.value].from_config(config, start)

    @property
    def nmodes(self):
        return self._db.nmodes

    @property
    def dynsampling(self):
        return self._db.dynsampling

    @property
    def info(self):
        return self._db.info

    @property
    def molecule(self):
        return self._db.molecule

    @property
    def model_info(self):
        return self._db._model_info

    @property
    def system(self):
        return self._db.system

    @property
    def increase(self):
        return self._db.increase

    @property
    def natoms(self):
        return self.molecule.natoms

    @property
    def atomids(self):
        return self.molecule.atomids

    @property
    def method(self):
        return self._db["method"]

    @property
    def modes(self):
        return self._db.modes

    @property
    def model(self):
        return self._db.model

    @property
    def nconditions(self):
        return len(self._db["crd"])

    @property
    def masses(self):
        return self._db.masses
