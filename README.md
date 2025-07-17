# Dipolar-Exchange Spin Waves in Magnetic Multilayers Incorporating Antimagnons

This repository provides **MATLAB** and **Python** implementations for calculating and visualising dipolar-exchange spin-wave spectra in a variety of magnetic multilayers, including the effect of antimagnons.

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

* Full explenation in Ref. [1] (`arXiv:2412.10888`) and benchmarks against the exchange–dipolar model of Ref. [2] (`Phys. Rev. Lett. 128, 217201`)
* Works in  MATLAB 
* Arbitrary number of layers, each with user-defined  
  * thickness  
  * saturation magnetisation  
  * exchange stiffness  
  * gyromagnetic ratio  
* Includes long-range dipolar coupling **and** intra-layer exchange automatically.  
* On-the-fly sweeping of effective magnetic field **H** overcomed by each layer.  
* Generates publication-quality band-structure plots.

---

## Repository Structure

| Path | Contents |
|------|----------|
| `matlab/` | `*.m` scripts and functions |
| `examples/` | Ready-to-run demo scripts for each system |
| `docs/` | Auto-generated API reference (Sphinx / MATLAB help) |
| `data/` | Sample material parameter files (`.json`, `.mat`) |

---

## Requirements

### MATLAB
* R2020b (or newer)  
* `Signal Processing Toolbox` (for plotting utilities)

