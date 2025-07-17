# Dipolar-Exchange Spin Waves in Magnetic Multilayers Incorporating Antimagnons

This repository provides **MATLAB** and **Python** implementations for calculating and visualising dipolar-exchange spin-wave spectra in a variety of magnetic multilayers, including the effect of non-equilibrium antimagnons.

---

## Implemented Modules

1. **Band-structure solvers**  
   a. Antiparallel ferromagnetic multilayer  
   b. Parallel ferromagnetic multilayer  
   c. Magnetostatic surface spin waves (MSSWs) in a semi-infinite ferromagnetic bulk  
   d. Antiferromagnetic / ferromagnetic (AFM/FM) hybrid multilayer  

2. **Tunability demo**  
   • Coherent/dissipative inter-layer coupling and layer-resolved chirality in an antiparallel FM bilayer.

---

## Features

* Fully reproduces the theory of Ref. [1] (`arXiv:2412.10888`) and benchmarks against the exchange–dipolar model of Ref. [2].  
* Works in both MATLAB and Python (NumPy/SciPy) with near-identical APIs.  
* Arbitrary number of layers, each with user-defined  
  * thickness  
  * saturation magnetisation  
  * exchange stiffness  
  * Gilbert damping  
  * gyromagnetic ratio  
* Includes long-range dipolar coupling **and** intra-layer exchange automatically.  
* On-the-fly sweeping of external magnetic field **H** and spin-orbit torque **τ** vectors.  
* Generates publication-quality band-structure plots (`.png`, `.pdf`).

---

## Repository Structure

| Path | Contents |
|------|----------|
| `matlab/` | `*.m` scripts and functions |
| `python/` | `*.py` modules and Jupyter notebooks |
| `examples/` | Ready-to-run demo scripts for each system |
| `docs/` | Auto-generated API reference (Sphinx / MATLAB help) |
| `data/` | Sample material parameter files (`.json`, `.mat`) |

---

## Requirements

### MATLAB
* R2020b (or newer)  
* `Signal Processing Toolbox` (for plotting utilities)

### Python
* Python ≥ 3.8  
* `numpy`, `scipy`, `matplotlib`, `numba` (optional speed-up)  
* `jupyter` (for the tutorial notebooks)

Install with:

```bash
pip install -r python/requirements.txt
