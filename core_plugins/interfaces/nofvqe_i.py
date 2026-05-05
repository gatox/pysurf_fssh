import numpy as np

from nofvqe import NOFVQE
from pysurf.spp.qm import AbinitioBase
from jinja2 import Template #build templates
from pysurf.system import Molecule

#geometry in bohr units
tpl = Template("""
unit bohr
{{chg}} {{mult}} {% for atomid, crd in mol %} 
 {{mol.format(atomid, crd)}} {% endfor %}

""")

class IntNOFVQE(AbinitioBase):
    """Interface for the quantum algorithm NOFVQE.
    """

    _user_input = """
    #==========================================================================
    #                     Interface for NOFVQE with PySurf
    #==========================================================================
    #--------------------------------------------------------------------------
    # Charge:
    #--------------------------------------------------------------------------
    chg = 0 :: int
    #--------------------------------------------------------------------------
    # Multiplicity:
    #--------------------------------------------------------------------------
    mult = 1 :: int
    #--------------------------------------------------------------------------
    # Basis set:
    #--------------------------------------------------------------------------
    basis = sto-3g :: str :: sto-3g, 6-31G, cc-pVDZ
    #--------------------------------------------------------------------------
    # Functional:
    #--------------------------------------------------------------------------
    functional = pnof4 :: str :: pnof4, pnof5, pnof7, pnof8, vqe
    #--------------------------------------------------------------------------
    # Only Double pair:
    #--------------------------------------------------------------------------
    pair_double = False :: bool 
    #--------------------------------------------------------------------------
    # Convergent tolerance:
    #--------------------------------------------------------------------------
    conv_tol = 1.0e-6 :: float 
    #--------------------------------------------------------------------------
    # Max iterations:
    #--------------------------------------------------------------------------
    max_iterations = 500 :: int  
    #--------------------------------------------------------------------------
    # Optimization circuit methods:
    #--------------------------------------------------------------------------
    opt_circ = sgd :: str :: sgd, adam, spsa, cobyla, slsqp, l-bfgs-b, cmaes
    #--------------------------------------------------------------------------
    # Gradient: Ground state only available
    #--------------------------------------------------------------------------
    gradient = :: str :: analytics, df_fedorov, df_normal  
    [gradient(analytics)]
    analytics_grad = True :: bool
    #--------------------------------------------------------------------------
    [gradient(df_fedorov)]
    #--------------------------------------------------------------------------
    # Displaced geometries:
    #--------------------------------------------------------------------------
    d_shift = 1.0e-3 :: float 
    [gradient(df_normal)]
    #--------------------------------------------------------------------------
    # Displaced geometries:
    #--------------------------------------------------------------------------
    d_shift = 1.0e-3 :: float
    #--------------------------------------------------------------------------
    # Device: Note that the noise_simulator and real options only work with an
    #         IBM account. The hybrid options use the simulator. When 
    #         the optimal value is reached, the real QC or noise_simulator
    #         recomputes it.
    #--------------------------------------------------------------------------
    device = :: str :: simulator, noise_simulator, real, hybrid_real, hybrid_noise_simulator
    [device(simulator)]
    dev_simulator = True :: bool
    #--------------------------------------------------------------------------
    [device(hybrid_real)]
    #--------------------------------------------------------------------------
    # Number of shots: Circuit measurements for expectation values.
    # Optimization level: Degree of circuit transpilation.
    # Resilience level: Degree of error mitigation.
    #--------------------------------------------------------------------------
    n_shots = 1000 :: int 
    optimization_level = 0 :: int
    resilience_level = 0 :: int
    #--------------------------------------------------------------------------
    [device(hybrid_noise_simulator)]
    #--------------------------------------------------------------------------
    # Number of shots: Circuit measurements for expectation values.
    # Optimization level: Degree of circuit transpilation.
    # Resilience level: Degree of error mitigation.
    #--------------------------------------------------------------------------
    n_shots = 1000 :: int 
    optimization_level = 0 :: int
    resilience_level = 0 :: int
    #--------------------------------------------------------------------------
    [device(noise_simulator)]
    #--------------------------------------------------------------------------
    # Number of shots: Circuit measurements for expectation values.
    # Optimization level: Degree of circuit transpilation.
    # Resilience level: Degree of error mitigation.
    #--------------------------------------------------------------------------
    n_shots = 1000 :: int 
    optimization_level = 0 :: int
    resilience_level = 0 :: int
    #--------------------------------------------------------------------------
    [device(real)]
    #--------------------------------------------------------------------------
    # Number of shots: Circuit measurements for expectation values.
    # Optimization level: Degree of circuit transpilation.
    # Resilience level: Degree of error mitigation.
    #--------------------------------------------------------------------------
    n_shots = 1000 :: int
    optimization_level = 0 :: int
    resilience_level = 0 :: int
    """
    tpl = tpl

    implemented = ['energy', 'gradient','parameter','rdm1_opt', 'n_opt', 'vecs_opt']

    def __init__(self, 
                 config, 
                 atomids, 
                 nstates, 
                 basis, 
                 chg, 
                 mult,
                 functional,
                 conv_tol,
                 max_iterations,
                 opt_circ,
                 gradient,
                 pair_double,
                 d_shift,
                 device,
                 n_shots,
                 optimization_level,
                 resilience_level,
                 ):
        self.molecule = Molecule(atomids, None)
        self.natoms = len(atomids) 
        self.nstates = nstates
        self.basis = basis
        self.chg = chg
        self.mult = mult
        self.functional = functional
        self.conv_tol = conv_tol
        self.init_param = None
        self.max_iterations = max_iterations
        self.opt_circ = opt_circ
        self.gradient = gradient
        self.pair_double=pair_double
        self.C_MO = None
        if self.gradient == "analytics":
            self.d_shift = None
        else:
            self.d_shift = config["gradient"]["d_shift"]
        self.device = device
        if self.device == "simulator":
            self.n_shots = None
            self.optimization_level = None
            self.resilience_level = None
        else:
            self.n_shots = config["device"]["n_shots"]
            self.optimization_level = config["device"]["optimization_level"]
            self.resilience_level = config["device"]["resilience_level"]
        self.count = 1
        self.count_2 = 1

    @classmethod
    def from_config(cls, 
                    config, 
                    atomids, 
                    natoms, 
                    nstates=None, 
                    **kwargs):
        return cls(config, 
                   atomids,
                   nstates, 
                   config['basis'], 
                   config['chg'], 
                   config['mult'], 
                   config['functional'], 
                   config['conv_tol'],
                   config['max_iterations'],
                   config['opt_circ'],
                   config['gradient'],
                   config['pair_double'],
                   config.get("d_shift",None),
                   config['device'],
                   config.get("n_shots",None),
                   config.get("optimization_level",None),
                   config.get("resilience_level",None),
                   )


    def get(self, request):
        
        # Update coordinates
        self.molecule.crd = request.crd

        self._do_nofvqe_ene_grad()


        self.count_2 +=1

        # Output requested properties
        if 'energy' in request:
            self._out_energy(request)
        if 'gradient' in request:
            self._out_gradient(request)
        if 'parameter' in request:
            self._out_parameter(request)
        if 'rdm1_opt' in request:
            self._out_rdm1(request)
        if 'n_opt' in request:
            self._out_n(request)
        if 'vecs_opt' in request:
            self._out_vecs(request)
            self.count +=1
        return request

    def _do_nofvqe_ene_grad(self):
        print(f"nofvqe_i.py Get function is called {self.count_2} times")
        print("nofvqe_i.py Crd:", self.molecule.crd)
        string_geo = self.tpl.render(chg=self.chg, mult=self.mult,
                     mol=self.molecule)
            
        nofvqe_class = NOFVQE(string_geo,
                      functional=self.functional,
                      conv_tol=self.conv_tol,
                      init_param=self.init_param,
                      basis=self.basis,
                      max_iterations=self.max_iterations,
                      opt_circ=self.opt_circ,
                      gradient=self.gradient,
                      pair_double=self.pair_double,
                      d_shift=self.d_shift,
                      C_MO = self.C_MO,
                      dev=self.device,
                      n_shots=self.n_shots,
                      optimization_level=self.optimization_level,
                      resilience_level=self.resilience_level,
                      )
        E_min, params_opt, rdm1_opt, n_opt, vecs_opt, cj12, ck12, C_opt, elag = nofvqe_class.run_scnofvqe()
        # Saving optimal variables for the next iteration
        self.params_opt = params_opt
        self.init_param = params_opt
        self.C_MO = C_opt
        
        print("nofvqe_i.py After called _do_nofvqe_ene_grad and computed params_opt:",self.init_param)
        """Saving energy and gradient for the ground state"""
        self.energy = E_min
        self.rdm1_opt = rdm1_opt
        self.n_opt = n_opt
        self.vecs_opt = vecs_opt
        #self.grad = nofvqe_class.grad()
        self.grad = nofvqe_class._nuclear_gradient_analytics(n_opt,C_opt,cj12,ck12,elag)

    def _out_energy(self, request):
        print(f"nofvqe_i.py Function energy is called {self.count} times")
        """Energy of the ground state"""
        out_ene = self.energy
        request.set('energy', out_ene)

    def _out_parameter(self, request):
        print(f"nofvqe_i.py Function parameter is called {self.count} times")
        """Optimal parameter"""
        out_parameter = self.params_opt
        request.set('parameter', out_parameter)
    
    def _out_rdm1(self, request):
        print(f"nofvqe_i.py Function rmd1_opt is called {self.count} times")
        """Optimal rdm1"""
        out_rdm1 = self.rdm1_opt
        request.set('rdm1_opt', out_rdm1)

    def _out_n(self, request):
        print(f"nofvqe_i.py Function n_opt is called {self.count} times")
        """Optimal n"""
        out_n = self.n_opt
        request.set('n_opt', out_n)

    def _out_vecs(self, request):
        print(f"nofvqe_i.py Function vecs_opt is called {self.count} times")
        """Optimal vecs"""
        out_vecs = self.vecs_opt
        request.set('vecs_opt', out_vecs)

    def _out_gradient(self, request):
        print(f"nofvqe_i.py Function gradient is called {self.count} times")
        """Gradientof the ground state"""
        out_gradient = {}
        for state in request.states:
            if state == 0:
                out_gradient[state] = np.array(self.grad)
            else:
                raise SystemExit("Gradient for excited states have not yet been implemented") 
        request.set('gradient', out_gradient)

if __name__=='__main__':
    from pysurf.database import PySurfDB
    from pysurf.spp.request import Request
    from numpy import copy 

    db_file = "sampling.db"
    db = PySurfDB.load_database(db_file, read_only=True)
    crd = copy(db['crd'][0])
    atomids = copy(db['atomids'])
    natoms = len(crd)

    out = IntNOFVQE.from_questions(config="spp_main.inp",
                                   atomids=atomids,
                                   natoms=natoms,
                                   nstates=None)
    
    # Create request for ground state (state 0)
    request = Request(crd, ['energy', 'gradient'], [0])

    # Run calculation
    response = out.get(request)

    # Print results
    print("Energy:", response['energy'])
    print("Gradient:", response['gradient'][0])
