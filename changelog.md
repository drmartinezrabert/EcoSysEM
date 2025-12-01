# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased] - yyyy/mm/dd
### Added
- Bioenergetic calculations.
  - Catabolic cell-specific power.
  - Empirical requirement cell-specific power.
- New environment modelling framework.
  - Multi-Metabolic State Model (MMSM) based on cell-specific power.
### Changed
- README file.

## [0.4] - yyyy/mm/dd
### Added
- Local sensitivity analysis of Gibbs free energy: `ThSA.local_sa_DeltaGr()` and `Environment.local_sa_DGr()`.
  - Local (Derivative).
  - Sigma-normalized derivative.
  - Reference-normalized derivative.
  - Variance-normalized derivative.
  - Pearson's correlation. 
- New environment in `Hydrosphere`.
  - Water column: `WaterColumn()`.
- New arguments in plotting functions: `plotVarMap2D()`, `plotZonalMean()` and `plotCrossSections()`.
  - Font family: `fontFamily='Arial'`.
  - Font weight of title: `fwtl='normal'`
<!--
- Global sensitivity analysis of Gibbs free energy: `ThSA.global_sa_DeltaGr()` and `Environment.global_sa_DGr()`.
  - Sobol' indices.
-->
### Fixed
- Bug in `ThEq.ThEq.pHSpeciation()`: Now the function handle nan values when all caompounds are requested (`rAllConc = True`).
### Changed
- README file.
- Rename the '2D compound sensitivity analysis (contourf plot)'.
  - `Environment.conc_var_DGr()` and `ThSA.conc_var_DeltaGr()` instead of `Environment.conc_sa_DGr()` and `ThSA.conc_sa_DeltaGr()`.
- `Environment.getDGr()`. 'Hydrosphere.WaterColumn()' environment has been incorporated.
- Improvement of plotting - tick labels: `plotVarMap2D()`, `plotZonalMean()` and `plotCrossSections()`.
- Design improvment on `plotCrossSections()`.

## [0.3] - 2025/11/07
### Added
- New environment in `Hydrosphere`.
  - General (or non-specific) Water Body: `GWB()`.
- Thermodynamic State Analysis (ThSA) is incorporated in `Environment class`.
  - ThSA functions are behaviours of `ISA()`, `ISAMERRA2()`, `CAMSMERRA2()` and `GWB()`.
- Sensitivity analysis of non-standard Gibbs free energy.
  - 2D compound sensitivity analysis (contourf plot).
- Creation of plotting script
  - New function: plotVarMap2D()
  - New function: plotZonalMean()
  - New function: plotCrossSections()
### Fixed
- Bug in `Environment.loadData()`: data existence wasn't check before loading.
- Bug in `CAMSMERRA2()`: error came out when phase was defined.
### Changed
- README file.

## [0.2] - 2025/11/03
### Added
- New atmospheric models (environment definition).
  - Modern-Era Retrospective analysis for Research and Applications, Version 2 (MERRA-2).
  - ISA-MERRA2 atmospheric model.
  - Copernicus Atmosphere Monitoring Service (CAMS).
  - CAMS-MERRA2 atmospheric model.
- New functions to calculate kinetic rates.
  - Michaelis-Menten.
  - Michaelis-Menten + Arrhenius .
- Calculation of non-standard Gibbs free energy (non-ideal liquids).
  - Estimation of thermodynamic activity coefficients.
    - Activity coefficients of electrolytes. Debye-Hückel theory (extended version).
    - Activity coefficients of non-electrolyte in electrolyte solutions. Setschenow-Schumpe equation.
- Module to run **EcoSysEM platform** via Command Line Interface (CLI) to download data: getDataMERRA2() and getDataCAMS().
### Changed
- README file.
- Update thermodynamic, equilibrium and rate functions for 3D data.
  - pH speciation.
  - Solubility (Henry's law).
  - Gibbs free energy calculations.
  - Kinetics.

## [0.1] - 2024/10/16
First pre-release of EcoSysEM platform (v0.1). 
### Added
- Environment definition.
    - International Standard Atmosphere model (ISA).
- Thermodynamic state analysis (ThSA).
    - Calculation of non-standard Gibbs free energy (ideal fluid).
- Equilibrium calculations.
    - pH speciation.
    - Solubility (Henry's law).
    - Estimation of equilibrium constants (Keq).
- Framework based on `.csv` files to define reactions.
- Framework based on `.csv` files to define constants (Henry's solubility constants) and thermodynamic parameters (DG0f, DH0f...).
- Proper README file.
<!--
### Added
- Lorem ipsum...
### Changed
- Lorem ipsum...
### Fixed
- Lorem ipsum...
### Removed
- Lorem ipsum...
-->
