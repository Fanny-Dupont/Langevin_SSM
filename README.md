# Reproducible Code — Simulation & Case Study

This repository contains all reproducible code for the simulation study and 
case study presented in:

**Inferring resource selection and utilization distributions from irregular 
and error-prone animal tracking data**  
Fanny Dupont, Brett McClintock, Jan-Ole Fischer, Marianne Marcoux, 
Nigel E. Hussey, Marie Auger-Méthé

---

## Repository Structure

- `simulation_study/` — Code to reproduce the simulation study. Uses the [`langevinSSM`](https://github.com/bmcclintock/langevinSSM) package.
- `Code_Case_Study/` — Code to reproduce the case study. Uses both the [`langevinSSM`](https://github.com/bmcclintock/langevinSSM) package and hand-coded functions and Cpp template.
- `Additional_simulation_penalty/` — Additional simulation results, implemented using hand-coded functions.

---

## Package

Functions for fitting the Langevin SSM model are available in the 
[`langevinSSM`](https://github.com/bmcclintock/langevinSSM) package. 
See the package repository for installation instructions and documentation.

---

## Contact

For questions, please contact [Fanny Dupont](fanny.dupont@stat.ubc.ca) .
