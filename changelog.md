# Changelog
All notable changes to this project will be documented in this file.

<!--
## [Unreleased] - yyyy/mm/dd
### Added
- _Lorem ipsum..._
### Removed
- _Lorem ipsum..._
### Fixed
- _Lorem ipsum..._
### Changed
- _Lorem ipsum..._
-->

## [0.6] - yyyy/mm/dd
### Added
- New optional arguments in `WaterColumn.plotVariables()`: `colors`, `title`, `title_fs`.
- New optional arguments in `Environments.getDGr()`: `printDG0r`,`printDH0r`.
### Removed
- _Lorem ipsum..._
### Fixed
- Bug in `ThSA.ionicStrength()`. Now it handels NaN value in compound concentration. 
### Changed
- README file.

## [0.5.1] - 2026/03/16
### Added
- New functions in `thermodynamics.py`:
  - `ThSA.getDeltaHr()`.
  - `ThEq.get_concentrations_Henry()`.
- New function in `environments.py`: `Environment.getDHr()`.
### Fixed
- Bug in `Environment.getDGr()`: Argument 'specComp' wasn't optional. Now, default value is _False_.
- Bug in `ThSA.sobol_indices_DeltaGr`. A `TypeError` occurred when it shouldn't have.
### Changed
- README file.
- Update thermodynamic date.
  - Specific heat capacity.

## [0.5] - 2026/03/13
### Added
- Computation of contributions (weights) of each grid (defined by longitude, latitude and altitude arrays) to statistics (e.g., quantiles/quartiles). User can consider the weights with the argument `weights = True`, if statistics is involved in the function.
- New Command Line Interface (CLI) scripts:
  - `cmd_getDataMERRA2.py`.
  - `cmd_getDataCAMS.py`.
- New parser argument in `cmd_getDataCAMS.py`: `_variables`.
- New function for plotting variables of water column: `WaterColumn.plotVariables()`.
- New functions for estimating thermodynamic parameters of water: `ThP class` in `thermodynamics.py`.
  - Surface tension: `ThP.surface_tension()`.
  - Vapor pressure: `ThP.vapor_pressure()`.
  - Osmotic coefficient: `ThP.osmotic_coefficient()`.
  - Water activity: `ThP.water_activity()`.
- New function to plot seasonal variability of a variable in `plotting.py`: `plot_seasonality()`.
- New optional argument in `CAMS.getDataCAMS()` (also included in `cmd_getDataCAMS.py`): `drop_variables`.
- New optional arguments in `plotVarMap2D()`: `logColorbar`, `cb_minor_ticks`, `cb_ticks`, `num_cb_ticks`, `cb_label_rotation`, `projection`.
- New optional argument in `plotCrossSections()`: `logColorbar`, `cb_minor_ticks`, `cb_ticks`, `num_cb_ticks`, `cb_label_rotation`, `projection`.
- New optional argument in `CAMSMERRA2`: `fillMissing`.
- New optional arguments in `plotZonalMean()`: `cb_ticks`, `num_cb_ticks`, `numlevels`, `grid_weights`.
- New function for `CAMS()` model: `fill_missing_levels()`. Fill in the missing data for the target levels with NaN values.
### Removed
- General Command Line Interface script: `ecosysem_cmd.py`. Now each function has its own script (`cmd_functionName.py`).
### Fixed
- Bug in `CAMSMERRA2()`: Missing key parameters when calling `CAMSMERRA2._interpolateCAMS()` (`target_lats` and `target_lons`).
- Bug in `cmd_getDataCAMS.py`: Now argument `_d` accepts string ('All'), integers and list of integers.
- Bugs in `Environment.combData()`:
  - Multiple years of CAMS data (>2 years) can be combined.
  - Multiple years months of MERRA2 data can be defined (before only first and last).
  - Only accepts valid models ('MERRA2' and 'CAMS').
### Changed
- README file.
- Density function moved to `ThP class`: `ThP.density()`.
- Visual improvement in plotting.

## [0.4] - 2025/12/17
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
  - Font weight of title: `fwtl='normal'`.
- Create new databases:
  - Specific heat capacities of compounds: `data/Cpi.csv`.
  - Equilibrium equation of electrolites for activity estimations: `reactions/electrolytes.csv`.
- Computation of heat capacity of reaction from specific heat capacity of compounds (C<sub>pi</sub>): `ThP.getDeltaCp()`.
- Global sensitivity analysis of Gibbs free energy.
  - Variance-based sensitivity analysis or Sobol’ indices: `ThSA.sobol_indices_DeltaGr()` and `Environment.sobol_indices_DGr()` for `GWB()`.
- Bioenergetic calculations: `bioenergetics.py` and `Environment` behaviours.
  - Catabolic cell-specific power.
  - Empirical requirement cell-specific power.
- New environment modelling framework.
  - Multi-Metabolic State Model (MMSM) based on cell-specific power.
### Fixed
- Bug in `ThEq.pHSpeciation()`: Now the function handle nan values when all caompounds are requested (`rAllConc = True`).
- Bug in `ThP.getThP()`: Now the function return empty Param if parameter(s) are not available for a compound.
- Bug in `ThSA.getDeltaGr()`: Now the function handle HCO3- or other speciation compounds of CO2 as a reaction product.
- Warning in `ThSA.local_sa_DeltaGr()`. Deprecation NumPy 1.25 solved. More information [here](https://github.com/pybamm-team/PyBaMM/issues/3052).
### Changed
- README file.
- Rename the '2D compound sensitivity analysis (contourf plot)'.
  - `Environment.conc_var_DGr()` and `ThSA.conc_var_DeltaGr()` instead of `Environment.conc_sa_DGr()` and `ThSA.conc_sa_DeltaGr()`.
- `Environment.getDGr()`. 'Hydrosphere.WaterColumn()' environment has been incorporated.
- Improvement of plotting - tick labels: `plotVarMap2D()`, `plotZonalMean()` and `plotCrossSections()`.
- Design improvment on `plotCrossSections()`.
- Improvement of Gibbs free energy calculation: `ThSA.getDeltaGr()`.
  - Non-standard enthalpy of reaction in function of temperature and heat capacity of reactions is included. 

## [0.3] - 2025/11/07
### Added
- New environment in `Hydrosphere`.
  - General (or non-specific) Water Body: `GWB()`.
- Thermodynamic State Analysis (ThSA) is incorporated in `Environment class`.
  - ThSA functions are behaviours of `ISA()`, `ISAMERRA2()`, `CAMSMERRA2()` and `GWB()`.
- Sensitivity analysis of non-standard Gibbs free energy.
  - 2D compound sensitivity analysis (contourf plot).
- Creation of plotting script.
  - New function: plotVarMap2D().
  - New function: plotZonalMean().
  - New function: plotCrossSections().
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
  - Michaelis-Menten + Arrhenius.
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
