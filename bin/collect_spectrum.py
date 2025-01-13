import os

#from copy_execute import CopyExecute
from setup_spectrum import SetupSpectrum
from colt import Colt
from pysurf.utils import exists_and_isfile
from pysurf.utils import SubfolderHandle
from pysurf.logger import get_logger, Logger
from pysurf.sampling import Sampling


class CollectSpectrum(SubfolderHandle, Colt):
    """ Class to collect and evaluate spectral information from the single point calculations

        It uses folder and subfolder from SetupSpectrum
    """
    folder = SetupSpectrum.folder
    subfolder = SetupSpectrum.subfolder
    file = 'init.db'
    props = ['energy', 'fosc', 'crd']

    _user_input = """
    sampling_input = sampling.inp :: existing_file
    """

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        self.logger = get_logger('collect_spectrum.log', 'collect_spectrum.log')
        subfolderhandle = SubfolderHandle(self.folder, self.subfolder)
        specfile = os.path.join(self.folder, 'spectrum.db')

        if exists_and_isfile(specfile):
            self.sampling = Sampling.from_db(config['sampling_input'], specfile, logger=self.logger)
            self.dimensions = self.sampling.info['dimensions']
            self.variables = self.sampling.info['variables']
            if not self._check_sampling(self.sampling):
                logger.error(f"Existing spectrum db is corrupted")
            for counter in range(self.sampling.nconditions, len(subfolderhandle)):
                snew = Sampling.from_db(config['sampling_input'], subfolderhandle.get_file(self.file, counter), logger=self.logger, sp=True)
                if self._check_sampling(snew):
                    self.add_condition(snew)

        else:
            counter = 0
            for file in subfolderhandle.fileiter('init.db'):
                snew = Sampling.from_db(config['sampling_input'], file, logger=self.logger, sp=True)
                if counter == 0:
                    self.sampling = Sampling.create_db(config['sampling_input'], specfile, snew.info['variables'], snew.info['dimensions'],
                                                    snew.molecule, snew.modes, snew.model, sp=False, logger=self.logger, fillit=False)
                    self.dimensions = self.sampling.info['dimensions']
                    self.variables = self.sampling.info['variables']
                self.add_condition(snew)
                counter += 1

    def add_condition(self, snew):
        for prop in self.props:
            self.sampling.append(prop, snew.get(prop, 0))
        self.sampling.increase

    def _check_sampling(self, samp):
        if not all(item in samp.info['variables'] for item in self.props):
            return False
        if samp.info['dimensions']['natoms'] != self.dimensions['natoms']:
            return False
        return True


if __name__ == "__main__":
    CollectSpectrum.from_commandline()
