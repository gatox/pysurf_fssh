import numpy as np

from ..database import PySurfDB
from colt import Colt
import os


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
    # Substeps: True or False
    # Nore: This feature is only available for Surface_Hopping. 
    #       For other methods, press Enter to continue.
    #--------------------------------------------------------------------------
    substeps = :: str, optional :: True, False 
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
    # Oscillation string (fosc), state-to-state transition moments (sts_mom)
    # and optimal parameter after vqe performance (parameter).
    # Note: The first two properties are available in Q-Chem and the latter 
    # one in NOFVQE.
    #       Press Enter to skip this step.
    #--------------------------------------------------------------------------
    save_properties = :: str, optional :: fosc, sts_mom, parameter
    #==========================================================================
    #                            Nose-Hoover thermostat
    #==========================================================================
    thermostat = :: str :: True, False
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
    rev_vel_no_hop = True :: bool :: True, False
    [substeps(True)]
    n_substeps = 10 :: int
    [substeps(False)]
    n_substeps = False :: bool
    [rescale_vel(momentum)]
    number_vdf = False :: str :: False, nonlinear, linear
    [rescale_vel(nacs)]
    res_nacs = True :: bool
    [thermostat(True)]
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
    [thermostat(False)]
    therm = False :: bool
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
            if config["substeps"] == "True":
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
        self.nob_dim = 0
        self.ekin = 0
        self.epot = 0
        self.grad = []
        self.nac = {}
        self.ene = []
        self.vk = []
        self.u = []
        self.rho = []
        
        if np.isscalar(self.mass):
            self.natoms = 1
        elif isinstance(self.mass, np.ndarray) != True:
            self.natoms = np.array([self.mass])
        if config["thermostat"] == "True":
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
            if config["save_properties"] == "parameter":
                self.save_properties += ["rdm1_opt", "n_opt", "vecs_opt"]
        else:
            self.save_properties = []

        self.additional = {}

    def save_additional(self, db):
        if not self.save_properties:
            return
        print("we are saving:", self.additional)
        for prop in self.save_properties:
            db.set(prop, self.additional[prop])


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

    @classmethod
    def from_db_frame(cls, db_file, config_file="prop.inp"):
        """
        Reconstruct a State from the last frame of a results.db file.
        Useful for restarting an interrupted trajectory.
        """
        if not os.path.exists(db_file):
            raise FileNotFoundError(f"Database file {db_file} not found.")

        # Load database and get last frame index
        db = PySurfDB.load_database(db_file, read_only=True)
        nframes = len(db["crd"])
        if nframes == 0:
            raise ValueError(f"No frames found in {db_file}")
        last = int(nframes - 1)

        # --- Reload configuration (prop.inp) ---
        state = cls.from_questions(config=config_file)  
        mdsteps_conf = int(state.mdsteps)
        if mdsteps_conf > last:
            # --- Reattach dynamic properties ---
            state.crd = np.copy(db["crd"][last])
            state.vel = np.copy(db["veloc"][last])
            state.grad = np.copy(db["gradient"][last]) if "gradient" in db else []
            state.ene = np.copy(db["energy"][last]) if "energy" in db else []
            state.t = float(np.copy(db["time"][last]))
            state.epot = float(np.copy(db["epot"][last]))
            state.ekin = float(np.copy(db["ekin"][last]))
            if state.method == "Surface_Hopping":
                state.instate = int(np.copy(db["currstate"][last]))
                state.nac = np.copy(db["nacs"][last]) if "nacs" in db else {}
                state.ncoeff = np.copy(db["populations"][last])

            print(f"[Restart] Loaded frame {last} from {db_file} at t = {state.t}")
            return state
        else:
            raise SystemExit(f"Last md iteration from results.db is {last} and md_steps from prop.inp is {mdsteps_conf}, so no need to restart :)")

if __name__ == "__main__":
    State.from_questions(config="prop.inp")
