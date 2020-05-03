import numpy as np

from pysurf.database import PySurfDB

class DynDB(PySurfDB):
    variables = ['crd_equi', 'modes_equi', 'model', 'atomids', 'freqs_equi', 'masses', 'currstate', 'crd', 'veloc', 'energy', 'ekin', 'epot', 'etot', 'time']

    @classmethod
    def from_dynamics(cls, dbfile):
        info = cls.info_database(dbfile)
        return cls.load_database(dbfile, info['variables'], info['dimensions'])

    @classmethod
    def create_db(cls, dbfile, sampling, nstates, props):
        variables = cls.variables
        for prop in props:
            if prop not in variables:
                variables += [prop]
        dims = sampling.info['dimensions']
        dims['nstates'] = nstates
        db = cls.generate_database(dbfile, variables, dims, model=sampling.model, sp=False)
        db.add_reference_entry(sampling.molecule, sampling.modes, sampling.model)
        return db

    def add_step(self, time, data, veloc, currstate, ekin, epot, etot):
        self.append('time', time)
        for entry in data:
            if entry == 'gradient':
                grad = np.empty((self.nstates, self.natoms, 3))
                for state in data[entry]:
                    grad[state, :, :] = data[entry][state]
                self.append(entry, grad)
                continue
            if entry in ['crd', 'fosc', 'energy', 'transmom']:
                self.append(entry, data[entry])
        self.append('currstate', currstate)
        self.append('ekin', ekin)
        self.append('epot', epot)
        self.append('etot', etot)
        self.increase


