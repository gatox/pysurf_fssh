# PySurf-FSSH: New Plugin for PySurf to Perform Nonadiabatic Dynamics using FSSH

# <img src="https://raw.githubusercontent.com/gatox/pysurf_fssh/master/docs/logo_pysurf_fssh.jpeg">

#

[PySurf-FSSH](https://doi.org/10.1021/acs.jctc.4c00012) is an extension of the [PySurf](https://github.com/mfsjmenger/pysurf) software package that enables nonadiabatic molecular dynamics simulations using the Fewest Switches Surface Hopping (FSSH) method. It includes built-in interfaces with several quantum chemistry packages such as [Q-Chem](https://manual.q-chem.com/latest/), [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas), [BAGEL](https://nubakery.org/), and the quantum algorithm State-Averaged Orbital-Optimized VQE ([SA-OO-VQE](https://github.com/qc2nl/SAOOVQE)). Thanks to PySurf’s flexible architecture, new interfaces can be added with minimal effort. As a demonstration, an analytical Linear Vibronic Coupling (LVC) model is also included, providing a fast and accessible way to test and extend PySurf-FSSH functionality.

---

# \:rocket: Installation

We recommend installing inside a fresh [Anaconda](https://www.anaconda.com/) environment:

```bash
conda create -y -n pysurf_saoovqe python=3.9
conda activate pysurf_saoovqe
```

Install dependencies for SAOOVQE:

```bash
pip install numpy scipy numba openfermion cirq
conda install -c psi4 psi4
```

Clone and install the SAOOVQE package:

```bash
git clone https://github.com/qc2nl/SAOOVQE.git
cd SAOOVQE
pip install -e .
cd ..
```

Clone and set up PySurf-FSSH:

```bash
git clone https://github.com/gatox/pysurf_fssh.git
cd pysurf_fssh
pip install pycolt qctools netcdf4 numpy matplotlib
```

Set environment variables (edit your `.bashrc` or `.zshrc`):

```bash
export PYTHONPATH=/full/path/to/pysurf_fssh:$PYTHONPATH
export SAOOVQEDIR=/full/path/to/SAOOVQE
```

---

# \:books: Purpose and Features

- Surface hopping dynamics via the Fewest Switches Surface Hopping (FSSH) method
- Support for ab initio on-the-fly dynamics (via Q-Chem, OpenMolcas, etc.)
- Support for model Hamiltonians (e.g., pyrazine LVC)
- Decoherence corrections (e.g., EDC)
- Velocity rescaling strategies (e.g., NACS-based)
- Thermostatting options for temperature control
- NetCDF-based trajectory database (`results.db`)

---

# \:sunrise: Examples

The `examples/` folder includes ready-to-run demonstrations:

## 1. CH2NH + SAOOVQE (Ab Initio)

Located in: `examples/ch2nh_saoovqe/`

Each subfolder (e.g., `ch2nh_f_noise_t_therm`) represents different combinations of noise and thermostat usage. Folder naming convention:

- `f` = false (not used)
- `t` = true (used)

For instance, `ch2nh_f_noise_t_therm` uses a thermostat, but no noise.

### Input setup:

- 10 trajectories, starting from the first excited state
- Time step \~0.25 fs
- Dynamics span 3 steps
- Temperature: 300 K
- Electronic structure settings: `spp.inp`
- Dynamics control: `prop.inp`

### How to run:

```bash
python ../../../pysurf_fssh/bin/sampling.py
python ../../../pysurf_fssh/bin/setup_propagation.py
python ../../../pysurf_fssh/bin/run_trajectory.py
```

Results are saved in `results.db` and `gen_results.out`.

## 2. Pyrazine LVC Model

Located in: `examples/pyrazine_lvc_model/`

Two folders:

- `pyrazine_f_therm`: no thermostat
- `pyrazine_t_therm`: thermostat enabled

### Input setup:

- 10 trajectories with Wigner sampling (`sampling.inp`)
- Time step: 20 a.u. (\~0.5 fs)
- Total time: 417 steps (\~200 fs)
- Settings in `spp.inp`, `prop.inp`

### How to run:

```bash
python ../../../pysurf_fssh/bin/sampling.py
python ../../../pysurf_fssh/bin/setup_propagation.py
python ../../../pysurf_fssh/bin/run_trajectory.py
```

This model runs extremely fast. To inspect results:

```bash
ncdump -v currstate results.db
cat gen_results.out
```

---

# \:hammer\_and\_wrench: Developer Notes

Clone repository:

```bash
git clone https://github.com/gatox/pysurf_fssh.git
cd pysurf_fssh
```

Optional: Save your current conda environment

```bash
conda env export > pysurf_saoovqe_env.yml
```

---

# \:card_index: Credits

This work was supported by the Innovational Research Incentives Scheme Vidi 2017 with project number **016.Vidi.189.044**, (partly) funded by the **Dutch Research Council (NWO)**.

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [`audreyr/cookiecutter-pypackage`](https://github.com/audreyr/cookiecutter-pypackage) project template.
