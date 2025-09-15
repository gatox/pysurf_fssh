# pysurf/sampling/fixed_sampler.py
import numpy as np
from colt import Colt

from ..system import Molecule
from ..system.atominfo import ATOMNAME_TO_ID, MASSES
from ..constants import U_TO_AMU

# SAMPLER base classes
from .base_sampler import DynSamplerBase, DynCondition


class Geometry(Colt):
    """
    Colt sub-block to provide geometry (either xyzfile OR atoms+coords).
    NOTE: Colt types:
      - atoms  -> list   (list of strings)
      - coords -> flist  (flat list of floats) or you can provide nested lists
      - velocities -> flist (flat list of floats) optional
    """
    _user_input = """
    xyzfile = :: str, optional
    atoms = :: list, optional
    coords = :: flist, optional
    velocities = :: flist, optional
    """

    @staticmethod
    def _read_xyz(path):
        atoms = []
        crd = []
        with open(path) as f:
            lines = [ln.rstrip() for ln in f if ln.strip()]
            # typical xyz: first line natoms, second comment, rest: atom x y z
            for line in lines[2:]:
                parts = line.split()
                atoms.append(parts[0])
                crd.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return atoms, np.array(crd, dtype=float)

    @classmethod
    def parse_block(cls, block):
        """
        Parse the Colt 'from' block (the dict that contains xyzfile / atoms / coords / velocities)
        Returns: (Molecule, velocities_array_or_None)
        """
        if block is None:
            raise ValueError("Geometry block is empty")

        # Case A: xyzfile provided
        if "xyzfile" in block and block["xyzfile"] is not None:
            atoms, crd = cls._read_xyz(block["xyzfile"])

        # Case B: inline atoms + coords
        elif "atoms" in block and block["atoms"] is not None and "coords" in block and block["coords"] is not None:
            atoms = block["atoms"]
            coords_val = block["coords"]

            # coords_val could be a flat list (flist) OR nested list (list of lists)
            crd = np.array(coords_val, dtype=float)
            if crd.ndim == 1:
                # flat list -> reshape
                if len(crd) % 3 != 0:
                    raise ValueError("coords length not divisible by 3")
                crd = crd.reshape(len(atoms), 3)
            elif crd.ndim == 2:
                # ensure shape matches atom count
                if crd.shape[0] != len(atoms) or crd.shape[1] != 3:
                    raise ValueError("coords shape doesn't match number of atoms")
            else:
                raise ValueError("coords array has unexpected number of dimensions")
        else:
            raise ValueError("Geometry requires either 'xyzfile' OR both 'atoms' and 'coords'")

        # Convert atom names -> atomic numbers (ids)
        try:
            atomids = np.array([ATOMNAME_TO_ID[a] for a in atoms], dtype=int)
        except Exception as e:
            raise ValueError(f"Error converting atom names to ids: {e}")

        # masses in a.u.
        masses = np.array([MASSES[idx] * U_TO_AMU for idx in atomids], dtype=float)

        # optional velocities
        velocities = None
        if "velocities" in block and block["velocities"] is not None:
            v = np.array(block["velocities"], dtype=float)
            if v.ndim == 1:
                if v.size != 3 * len(atoms):
                    raise ValueError("velocities length does not match 3*Natoms")
                v = v.reshape(len(atoms), 3)
            elif v.ndim == 2:
                if v.shape != (len(atoms), 3):
                    raise ValueError("velocities shape must be (N,3)")
            else:
                raise ValueError("velocities must be flat list or list of lists")
            velocities = v

        # Build Molecule object expected by PySurf (atomids, crd, masses)
        mol = Molecule(atomids, crd, masses)
        return mol, velocities


class FixedSampler(DynSamplerBase):
    """
    Deterministic sampler: returns a single condition with the provided geometry
    and zero (or user-provided) velocities.

    Usage in sampling.inp:

    [Sampling]
    method = FixedSampler
    sampling_db = sampling.db

    [FixedSampler]
    from = Geometry

    [Geometry]
    xyzfile = h2.xyz

    OR

    [Geometry]
    atoms = H H
    coords = 0.0 0.0 0.0  0.0 0.0 1.4
    velocities = 0.0 0.0 0.0  0.0 0.0 0.0   # optional
    """
    _user_input = """
    from = Geometry
    """

    _from = {"Geometry": Geometry}

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases(
            "from", {name: method.colt_user_input for name, method in cls._from.items()}
        )

    def __init__(self, system, velocities=None):
        # system is a Molecule instance
        self.system = system
        self.generated = False
        self._velocities = None if velocities is None else np.array(velocities, dtype=float)
        self.generated_count = 0
        self.nmax = 1  # default, overridden later
        
    @classmethod
    def from_config(cls, config, start=0):
        """
        config is the Colt-generated block for the sampler;
        the chosen 'from' block is nested inside config["from"].
        """
        from_block = config.get("from", None)
        if from_block is None:
            raise ValueError("FixedSampler.from_config: missing 'from' block")

        # Parse geometry and optional velocities using Geometry helper
        mol, velocities = Geometry.parse_block(from_block)
        return cls(mol, velocities)

    @classmethod
    def from_db(cls, database):
        # when building from an existing database, PySurf already has the reference system
        return cls(None)

    def get_init(self):
        return {"system": self.system, "modes": None}

    def get_condition(self):
        if self.generated_count >= self.nmax:
            return None
        if self.system is None:
            raise RuntimeError("FixedSampler: no system defined")
        crd = np.copy(self.system.crd)
        if self._velocities is None:
            veloc = np.zeros_like(crd, dtype=float)
        else:
            veloc = np.copy(self._velocities)
        state = 0
        self.generated = True
        self.generated_count += 1
        return DynCondition(crd, veloc, state)

    def get_number_of_conditions(self, nmax):
        self.nmax = nmax
        return nmax 
