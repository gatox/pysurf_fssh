from abc import abstractmethod
import os
import numpy as np
import time

from ..spp import SurfacePointProvider
from ..utils import exists_and_isfile
from ..logger import get_logger

from qctools.converter import Converter, time_converter
from colt import Plugin

from .dyn_db import DynDB


class PropagatorFactory(Plugin):
    _plugins_storage = "_methods"
    _is_plugin_factory = True

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases(
            "method", {name: method.questions for name, method in cls._methods.items()}
        )


class PropagatorBase(PropagatorFactory):
    _user_input = "inherited"
    extend_user_input: "inherited"

    _register_plugin = False
    # Properties have to be defined for each Propagator
    properties = ["energy", "gradient"]

    @abstractmethod
    def _run(self, nsteps, db, *args, **kwargs):
        """propagation routine"""
        pass

    def run(self, nsteps, dt, *args, **kwargs):
        self.output_header()
        self.start_time = time.perf_counter()
        self._run(nsteps, dt, *args, **kwargs)

    def get_runtime(self):
        return time.perf_counter() - self.start_time

    def __init__(
        self,
        spp_inp,
        sampling,
        nstates,
        nghost_states,
        restart=True,
        logger=None,
        properties=None,
    ):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model
        """
        self.nstates = nstates
        self.start_time = time.perf_counter()
        self.sampling = sampling

        if logger is None:
            self.logger = get_logger("prop.log", "propagator")
            info = {"spp inputfile": spp_inp}
            self.logger.header("PROPAGATION", info)
        else:
            self.logger = logger

        if properties is not None:
            for prop in properties:
                if prop not in self.properties:
                    self.properties += [prop]

        self.init = sampling.get_condition(0)
        # setup SPP
        if sampling.model is False:
            self.spp = SurfacePointProvider.from_questions(
                self.properties,
                nstates,
                sampling.natoms,
                nghost_states=nghost_states,
                atomids=sampling.atomids,
                config=spp_inp,
            )
        else:
            self.spp = SurfacePointProvider.from_questions(
                self.properties,
                nstates,
                sampling.nmodes,
                nghost_states=nghost_states,
                config=spp_inp,
            )

        if exists_and_isfile("prop.db"):
            self.db = DynDB.from_dynamics("prop.db")
            if len(self.db["crd"]) > 0:
                self.restart = True
            else:
                self.restart = False
        else:
            if self.init is None:
                self.logger.error(
                    "No initial condition or restart file for the Propagator"
                )
            else:
                self.create_new_db()
            self.restart = False

        if restart is False:
            self.restart = False

        if self.restart is True:
            self.masses = np.array(self.db["masses"])
            mode_output = "a"
        else:
            self.masses = sampling.masses
            mode_output = "w"

        self.t_converter = time_converter.get_converter(tin="au", tout="fs")
        self.output = get_logger("prop.out", "propagator_output", mode=mode_output)

    def create_new_db(self):
        name = "prop.db"
        if exists_and_isfile(name):
            os.remove(name)
        self.db = DynDB.create_db(name, self.sampling, self.nstates, self.properties)

    def output_header(self):
        self.output.info("#" + ("=" * 101))
        self.output.info(
            f"#{'Step':^9}|{'Time':^12}|{'State':^7}|{'Energy':^47}|{'Gradient':^11}|{'Runtime':^9}|"
        )
        self.output.info(
            f"#{' ':9}|{' ':^12}|{' ':^7}|{'kin':^11}|{'pot':^11}|{'tot':^11}|{'diff':^11}|{'RMS':^11}|{' ':^9}|"
        )
        self.output.info(
            f"#{' ':9}|{'[fs] ':^12}|{' ':^7}|{'[au]':^47}|{'[au]':^11}|{'[sec]':^9}|"
        )
        self.output.info("#" + ("=" * 101))

    def output_step(self, step, time, state, ekin, epot, etot, dE, grad=None):
        if grad is None:
            self.output.info(
                f"{step:>10}{self.t_converter(time):>13.2f}{state:>8} {ekin:>11.6f} {epot:>11.6f} {etot:>11.6f} {dE:>11.6f}{' ':12}{round(self.get_runtime(), 1):>10}"
            )
        else:
            grad = np.array(grad).flatten()
            rmsd = np.sqrt(grad.dot(grad) / grad.size)
            self.output.info(
                f"{step:>10}{self.t_converter(time):>13.2f}{state:>8} {ekin:>11.6f} {epot:>11.6f} {etot:>11.6f} {dE:>11.6f}{rmsd:11.6f}{round(self.get_runtime(), 1):>10}"
            )
