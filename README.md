# ESA Ramses Spacecraft Orbital Survival & SRP Perturbation Dynamics at (99942) Apophis (2029 Perigee)

[![License: MIT](https://shields.io)](https://opensource.org)
[![ORCID](https://shields.io)](https://orcid.org)

## Overview
High-precision 6-body numerical simulation framework developed in Python to evaluate the orbital stability and tidal survival of an **ESA Ramses-class rendezvous spacecraft** accompanying near-Earth asteroid **(99942) Apophis** across its April 13, 2029 terrestrial close approach.

## Physical & Numerical Framework
* **Integrator:** Implicit 15th-order Gauss–Radau adaptive solver (`IAS15`, local truncation tolerance $\epsilon = 10^{-12}$) implemented within the `REBOUND` computational ecosystem.
* **Ephemerides:** Initialized with live geocentric and heliocentric state vectors retrieved via `Astroquery` from NASA JPL Horizons (epoch April 10, 2029).
* **Hill Sphere Collapse & Tidal Stripping:** Models the extreme contraction of Apophis's Earth-dominated Hill sphere radius ($r_H \approx 0.57\text{ km}$ at perigee), contrasting surviving tight station-keeping orbits ($r_0 = 1.2\text{ km}$) against tidally stripped outer orbits ($r_0 = 5.0\text{ km}$).
* **Solar Radiation Pressure (SRP) Critical Sweep:** Evaluates 20 initial station-keeping altitudes ($r_0 \in [0.6, 4.4]\text{ km}$) under continuous solar photon acceleration ($A/m = 0.01\text{ m}^2/\text{kg}$, $C_R = 1.3$) to determine the maximum safe survival radius during perigee passage.

## Repository Structure
* `ramses_apophis_perigee_survival.py`: Coupled 6-body baseline simulation comparing bound vs. tidally stripped spacecraft trajectories.
* `ramses_srp_critical_sweep.py`: Parametric altitude sweep incorporating non-gravitational Solar Radiation Pressure (SRP) forces on spacecraft solar arrays.

## Author
**René Sagal Andrade**  
Departamento de Astrodinámica y Ciencias Computacionales, Universidad Nacional Autónoma de México (UNAM)  
ORCID: [0009-0003-5574-5384](https://orcid.org)

