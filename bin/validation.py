"""
PySurf Module:
    Validation and Training of Interpolators

Provide infrastructure for the training of interpolators
and test them against a validation set
"""
import numpy as np

from pysurf.database import PySurfDB
from pysurf.spp.dbinter import InterpolatorFactory
from pysurf.spp.request import RequestGenerator
from pysurf.logger import get_logger
from pysurf.colt import Colt


class Validation(Colt):

    _questions = """
    [validate]
    db = 
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_block("training", Training.questions)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        self.inter = Training.from_config(config['training'])
        self.inter.validate(config['validate']['db'], config['training']['properties'])


class Training(Colt):

    _questions = """
        database = :: existing_file
        properties = :: list
        interpolator = RbfInterpolator :: str
        energy_only = False :: bool
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("interpolator",
                                 {name: interpolator.questions
                                  for name, interpolator in InterpolatorFactory.interpolator.items()})

    @classmethod
    def from_config(cls, config):
        return cls(config['database'], config['properties'],
                   config['interpolator'], energy_only=config['energy_only'])

    def __init__(self, filename, properties, interpolator, energy_only=False):
        self.properties = properties
        self.db = PySurfDB.load_database(filename, read_only=True)
        self.logger = get_logger('db.log', 'validation', [])
        self.interpolator = InterpolatorFactory.plugin_from_config(interpolator,
                                                                   self.db,
                                                                   properties,
                                                                   logger=self.logger,
                                                                   energy_only=energy_only)

        self.interpolator.train()
        self.nstates = self.db.get_dimension_size('nstates')

    def validate(self, filename, properties):
        db = PySurfDB.load_database(filename, read_only=True)
        reqgen = RequestGenerator(self.nstates, properties, use_db=True)

        norm = {prop: [] for prop in properties}
        is_not_trustworth = 0
        ndata = len(db)

        for i, crd in enumerate(db['crd']):
            request = reqgen.request(crd, properties)
            result, is_trustworthy = self.interpolator.get(request)
            #
            for prop in properties:
                norm[prop].append(np.copy(result[prop] - db[prop][i]))
            #
            if not is_trustworthy:
                is_not_trustworth += 1

        for name, value in norm.items(): 
            self.compute_errors(name, value, nele)
        #
        trust = (ndata - is_not_trustworth)/ndata
        #
        self.logger.info(f"is_trustworthy = {trust}, number not trustworth={is_not_trustworth}")

    def compute_errors(self, name, prop, nele):
        prop = np.array(prop)
        #
        mse = np.sum(prop)/nele
        mae = np.sum(np.absolute(prop))/nele
        rmsd = np.sqrt(np.sum(prop**2)/nele)
        #
        maxval = np.amax(prop)
        minval = np.amin(prop)
        self.logger.info(f"{name}:\n mse = {mse}\n mae = {mae}\n"
                         f" rmsd = {rmsd}\n maxval = {maxval}\n minval={minval}\n")


def eucl_norm(x, y):
    return np.linalg.norm(x - y)


if __name__ == '__main__':
    Validation.from_commandline()
