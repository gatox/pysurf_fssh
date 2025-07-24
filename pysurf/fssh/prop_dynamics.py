import numpy as np

from ..database import PySurfDB
from colt import Colt


class State(Colt):

    _user_input = """ 
    #==========================================================================
    #               Parameters for Initialising Nuclear Propagation
    #==========================================================================
    db_file = :: existing_file 
    t = 0.0 :: float
    dt = 1.0 :: float
    mdsteps = 40000 :: float
    #--------------------------------------------------------------------------
    # Substeps: true or false
    # Nore: This feature is only available for Surface_Hopping. 
    #       For other methods, press Enter to continue.
    #--------------------------------------------------------------------------
    substeps = :: bool, optional :: true, false 
    #--------------------------------------------------------------------------
    # Rescale velocity: momentum or nacs
    # Note: This feature is only available for Surface_Hopping. 
    #       If 'momentum' is selected, the rescaling 
    #       could be improved by enabling 'number_vdf'.
    #       For other methods, press Enter to continue.
    #--------------------------------------------------------------------------
    rescale_vel = :: str, optional :: momentum, nacs 
    #--------------------------------------------------------------------------
    # instate is the initial state: 0 = G.S, 1 = E_1, ...
    #--------------------------------------------------------------------------
    instate = 0 :: int
    #--------------------------------------------------------------------------
    # Save additional properties: 
    # Oscillation string (fosc) and state-to-state transition moments (sts_mom)
    # Note: These properties are available in Q-Chem.
    #       Press Enter to skip this step.
    #--------------------------------------------------------------------------
    save_properties = :: str, optional :: fosc, sts_mom
    #==========================================================================
    #                             Nuclear Propagator
    #==========================================================================
    method = Born_Oppenheimer :: str :: Born_Oppenheimer, Surface_Hopping 
    [method(Born_Oppenheimer)]
    #==========================================================================
    #                         Ab_inito Molecular Dynamics
    #==========================================================================
    activate_BO = 0 :: int 
    [method(Surface_Hopping)]
    #==========================================================================
    #                       Nonadiabatic Molecular Dynamics 
    #==========================================================================
    nstates = 2 :: int
    states = 0 1 :: ilist
    ncoeff = 0.0 1.0 :: flist
    prob = tully :: str :: tully, lz, lz_nacs    
    coupling = nacs :: str :: nacs, wf_overlap, non_coup, semi_coup 
    decoherence = EDC :: str :: EDC, IDC_A, IDC_S, No_DC
    rev_vel_no_hop = true :: bool :: true, false
    [substeps(true)]
    n_substeps = 10 :: int
    [substeps(false)]
    n_substeps = false :: bool
    [rescale_vel(momentum)]
    number_vdf = false :: str :: false, nonlinear, linear
    [rescale_vel(nacs)]
    res_nacs = true :: bool
    #==========================================================================
    #                            Nose-Hoover thermostat
    #==========================================================================
    thermostat = :: bool :: true, false
    [thermostat(true)]
    #--------------------------------------------------------------------------
    # Friction coefficient
    #--------------------------------------------------------------------------
    xi = 0.0 :: float
    #--------------------------------------------------------------------------
    # Target tempertaure in Kelvin
    #--------------------------------------------------------------------------
    T = 300 :: float    
    #--------------------------------------------------------------------------
    # degrees of freedom 
    #--------------------------------------------------------------------------
    dof = nonlinear :: str :: nonlinear, linear
    [thermostat(false)]
    therm = false :: bool
    #==========================================================================
    """

    def __init__(
        self,
        config,
        crd,
        vel,
        mass,
        atomids,
        model,
        t,
        dt,
        mdsteps,
        instate,
        method,
        nstates,
        states,
        ncoeff,
        prob,
        rescale_vel,
        rev_vel_no_hop,
        coupling,
        decoherence,
        substeps,
        thermostat,
    ):
        self.crd = crd
        self.natoms = len(crd)
        self.atomids = atomids
        self.vel = vel
        self.mass = mass
        if model == 1:
            self.model = True
        else:
            self.model = False
        self.t = t
        self.dt = dt
        self.mdsteps = mdsteps
        self.instate = instate
        self.method = method
        if config["method"] == "Surface_Hopping":
            self.method = "Surface_Hopping"
            self.nstates = config["method"]["nstates"]
            self.states = config["method"]["states"]
            self.ncoeff = config["method"]["ncoeff"]
            self.prob = config["method"]["prob"]
            self.rescale_vel = config["rescale_vel"].value
            if config["rescale_vel"] == "momentum":
                self.reduced_kene = config["rescale_vel"]["number_vdf"]
            self.coupling = config["method"]["coupling"]
            if config["rescale_vel"] == "nacs":
                if self.coupling in ("wf_overlap", "non_coup"):
                    raise SystemExit(
                        "Incompatible coupling and rescaling"
                    )
            self.rev_vel_no_hop = config["method"]["rev_vel_no_hop"]
            self.decoherence = config["method"]["decoherence"]
            if config["substeps"] == "true":
                self.substeps = True
                self.n_substeps = config["substeps"]["n_substeps"]
            else:
                self.substeps = False
        elif config["method"] == "Born_Oppenheimer":
            self.method = "Born_Oppenheimer"
            self.activate_BO = config["method"]["activate_BO"]
        self.e_curr = None
        self.e_prev_step = None
        self.e_two_prev_steps = None
        self.ekin = 0
        self.epot = 0
        self.nac = {}
        self.ene = []
        self.vk = []
        self.u = []
        self.rho = []
        if np.isscalar(self.mass):
            self.natoms = 1
        elif isinstance(self.mass, np.ndarray) != True:
            self.natoms = np.array([self.mass])
        if config["thermostat"] == "true":
            self.thermostat = True
            self.xi = config["thermostat"]["xi"]
            self.dof = config["thermostat"]["dof"]
            self.t_target = config["thermostat"]["T"] * 3.166811e-6
            if self.dof == "nonlinear":
                self.q_eff = (3 * self.natoms - 6) * self.t_target * (10 * self.dt) ** 2
            else:
                self.q_eff = (3 * self.natoms - 5) * self.t_target * (10 * self.dt) ** 2
        else:
            self.thermostat = False

        if config["save_properties"] is not None:
            self.save_properties = [config["save_properties"]]
        else:
            self.save_properties = []

        self.additional = {}

    def save_additional(self, db):
        # either add fosc or sts_mom
        print("we are saving: ", self.additional)
        for prop, value in self.additional.items():
            db.set(prop, value)

    @classmethod
    def from_config(cls, config):
        crd, vel, mass, atomids, model = cls.read_db(config["db_file"])
        t = config["t"]
        dt = config["dt"]
        mdsteps = config["mdsteps"]
        instate = config["instate"]
        method = config["method"]
        nstates = config.get("nstates",None)
        states = config.get("states",None)
        ncoeff = config.get("ncoeff",None)
        prob = config.get("prob",None)
        rescale_vel = config.get("rescale_vel",None)
        rev_vel_no_hop = config.get("rev_vel_no_hop",None)
        coupling = config.get("coupling",None)
        decoherence = config.get("decoherence",None)
        substeps = config.get("substeps",None)
        thermostat = config.get("thermostat",None)
        return cls(
            config,
            crd,
            vel,
            mass,
            atomids,
            model,
            t,
            dt,
            mdsteps,
            instate,
            method,
            nstates,
            states,
            ncoeff,
            prob,
            rescale_vel,
            rev_vel_no_hop,
            coupling,
            decoherence,
            substeps,
            thermostat,
        )

    @staticmethod
    def read_db(db_file):
        db = PySurfDB.load_database(db_file, read_only=True)
        crd = np.copy(db["crd"][0])
        vel = np.copy(db["veloc"][0])
        atomids = np.copy(db["atomids"])
        mass = np.copy(db["masses"])
        model = np.copy(db["model"])
        if model == 1:
            model = True
        else:
            model = False
        return crd, vel, mass, atomids, model

    @classmethod
    def from_initial(
        cls,
        config,
        crd,
        vel,
        mass,
        atomids,
        model,
        t,
        dt,
        mdsteps,
        instate,
        method,
        nstates,
        states,
        ncoeff,
        prob,
        rescale_vel,
        rev_vel_no_hop,
        coupling,
        decoherence,
        substeps,
        thermostat,
    ):
        return cls(
            config,
            crd,
            vel,
            mass,
            atomids,
            model,
            t,
            dt,
            mdsteps,
            instate,
            method,
            nstates,
            states,
            ncoeff,
            prob,
            rescale_vel,
            rev_vel_no_hop,
            coupling,
            decoherence,
            substeps,
            thermostat,
        )


if __name__ == "__main__":
    State.from_questions(config="prop.inp")
