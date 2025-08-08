# Dipolar-Exchange Spin Waves in Magnetic Multilayers Incorporating Antimagnons

This repository provides **MATLAB** implementations for calculating and visualising dipolar-exchange spin-wave spectra in a variety of magnetic multilayers, including the effect of antimagnons.

---

## Implemented Modules

1. **Band-structure solvers**  
   • Antiparallel ferromagnetic multilayer (Fig. 2, Fig. S1)
   
   • Parallel ferromagnetic multilayer (Fig. 3, Fig. S2)
   
   • Magnetostatic surface spin waves (MSSWs) in a semi-infinite ferromagnetic bulk  (Fig. S3)
   
   • Antiferromagnetic / ferromagnetic (AFM/FM) hybrid multilayer  (Fig. 6)
   

3. **Tunability demo**  
   • Coherent/dissipative inter-layer coupling and layer-resolved chirality in an antiparallel FM bilayer. (Fig. 4) 

---

## Features

* Full explenation in Ref. [1] (`arXiv:2412.10888`) and aditional reference to the exchange–dipolar model of Ref. [2] (`Phys. Rev. Lett. 128, 217201`)
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

## Requirements

### MATLAB
* R2021b (or newer)  

