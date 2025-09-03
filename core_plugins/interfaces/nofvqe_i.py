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
    basis = sto-3g :: str :: cc-pvdz, sto-3g
    #--------------------------------------------------------------------------
    # Functional:
    #--------------------------------------------------------------------------
    functional = pnof4 :: str :: ca, ml, gu, bbc2, bbac3, cga, pnof4
    #--------------------------------------------------------------------------
    # Convergent tolerance:
    #--------------------------------------------------------------------------
    conv_tol = 1.0e-6 :: float 
    #--------------------------------------------------------------------------
    # Max iterations:
    #--------------------------------------------------------------------------
    max_iterations = 500 :: int  
    #--------------------------------------------------------------------------
    # Gradient: Ground state only available
    #--------------------------------------------------------------------------
    gradient = df_fedorov :: str :: df_fedorov, df_normal
    #--------------------------------------------------------------------------
    # Displaced geometries:
    #--------------------------------------------------------------------------
    d_shift = 1.0e-3 :: float 
    """
    tpl = tpl

    implemented = ['energy', 'gradient','parameter']

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
                 gradient,
                 d_shift):
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
        self.gradient = gradient
        self.d_shift = d_shift
        # self.icall = 0
        self._last_crd = None

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
                   config['gradient'],
                   config['d_shift'])


    def get(self, request):
        # if self.icall == 0:
        #     self.read_param = False
        #     self.icall = 1
        # else:
        #     self.read_param = True

        # Update coordinates
        self.molecule.crd = request.crd

        # Check if coordinates are the same as last call
        if self._last_crd is None or not np.allclose(self._last_crd, request.crd):
            self._do_nofvqe_ene_grad()
            self._last_crd = np.copy(request.crd)

        # Output requested properties
        if 'energy' in request:
            self._out_energy(request)
        if 'gradient' in request:
            self._out_gradient(request)
        if 'parameter' in request:
            self._out_parameter(request)
        return request

    def _do_nofvqe_ene_grad(self):
        string_geo = self.tpl.render(chg=self.chg, mult=self.mult,
                     mol=self.molecule)
                     
        nofvqe_class = NOFVQE(string_geo,
                      functional=self.functional,
                      conv_tol=self.conv_tol,
                      init_param=self.init_param,
                      basis=self.basis,
                      max_iterations=self.max_iterations,
                      gradient=self.gradient,
                      d_shift=self.d_shift)
        
        E_min, params_opt, _ = nofvqe_class.ene_vqe()
        self.init_param = params_opt
        """Saving energy and gradient for the ground state"""
        self.energy = E_min
        self.grad = nofvqe_class.grad()

    def _out_energy(self, request):
        """Energy of the ground state"""
        out_ene = self.energy
        request.set('energy', out_ene)

    def _out_parameter(self, request):
        """Energy of the ground state"""
        out_parameter = self.init_param
        request.set('parameter', out_parameter)

    def _out_gradient(self, request):
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

    out = IntNOFVQE.from_questions(config="spp.inp",
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
