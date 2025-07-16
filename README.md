# Binary Star System Simulator

A Python-based simulation and analysis of binary star systems, developed during research at Bradley University.


## Project Overview

This is an astronomy project that I started on May 2022 under the mentorship of Dr. Gregory Simonian at Bradley University. This project explores the orbital properties of rapidly rotating stars observed by the Kepler Telescope, testing the hypothesis that they are spectroscopic binary systems undergoing tidal locking. Developed during PHY 399 research at Bradley University under Dr. Gregory Simonian’s guidance, the simulation combines Newtonian mechanics with real Kepler data to generate and analyze radial velocity curves using Kepler's third law of planetary motion.

## Research Introduction

Inspired by a Kepler catalog of 27 unexpectedly fast-rotating stars, our study probes whether these objects form binary pairs bound by mutual gravity and tidal interactions. Spectroscopic binaries reveal stellar masses via Doppler shifts; by modeling radial velocity signatures, we can accept or reject our binary hypothesis.

## Research Methodology

- **Literature Review & Tools**: Studied spectroscopic binaries and Kepler’s laws; read Python Data Science Handbook. Used NumPy for mathematical functions, SciPy for integrators, Matplotlib for plotting, and AstroPy for astronomical constants, units, and data tables.
- **Model Development**: Derived amplitude and radial velocity equations via Kepler’s third law.
- **Data Integration**: Imported rotational periods (McQuillan et al., 2014) and masses (Berger et al., 2020).
- **Parameter Modeling**: Modeled inclination, mass ratio, and phase angle based on their statistical distributions; added white noise to account for observational errors.
- **Statistical Testing**: Generated thousands of synthetic radial velocity curves and used a Chi-square statistical test to identify best-fit binary parameters. We then compare our modeled radial velocity curve with the actual measured radial velocity times, measured by the Modspec instrument at the Michigan-Dartmouth-MIT (MDM) observatory.



## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Authors
- [@gregosimo](https://github.com/gregosimo)
- [@Christian0115](https://github.com/Christian-0115)
