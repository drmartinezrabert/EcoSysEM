# EcoSysEM platform
Eco-System Evaluation &amp; Modelling

*· Contributors: Eloi Martinez-Rabert, Begüm Nisa Kasaplı, Zoé Lemaire*.

*· Funded by:*

<table>
<tr><td><b>2024 - 2026:</b></td><td><img align="centre" src="https://github.com/user-attachments/assets/d03a5160-9186-46e3-8bf4-b55681ee1841" width="150"></td><td>RGY0058/2022</td></tr>
</table>

____________________________

## README Contents
- EcoSysEM platform | [GO](#ecosysem-platform)
- Before having fun... | [GO](#before-having-fun)
    - Anaconda Python installation | [GO](#gear-anaconda-python-installation)
    - Anaconda Navigator | [GO](#anaconda-navigator)
    - Anaconda Prompt or Terminal | [GO](#anaconda-prompt-or-terminal)
    - Spyder | [GO](#spyder)
    - Python packages | [GO](#python-packages)
    - Installation of packages using Anaconda Navigator | [GO](#installation-of-packages-using-anaconda-navigator)
    - Installation of packages using pip | [GO](#installation-of-packages-using-pip)
- Instructions for downloading and setting up EcoSysEM platform | [GO](#clipboard-instructions-for-downloading-and-setting-up-ecosysem-platform)
    - How to modify parameter databases | [GO](#how-to-modify-parameter-databases)
        - Database references | [GO](#database-references) 
    - How to modify reaction databases | [GO](#how-to-modify-reaction-databases)
    - Formulization of compounds | [GO](#formulization-of-compounds)
- Instructions to use EcoSysEM platform via Spyder | [GO](#clipboard-instructions-to-use-ecosysem-platform-via-spyder)
- EcoSysEM user guide | [GO](#ecosysem-user-guide)
    - EcoSysEM package layout | [GO](#ecosysem-package-layout)
    - Fundamentals and usage | [GO](#fundamentals-and-usage)
        - Environment definition and instance calling | [GO](#environment-definition-and-instance-calling)
        - Thermodynamic State Analysis (ThSA) | [GO](#thermodynamic-state-analysis-thsa)
        - Ecosystem modelling 🚧 | [GO](#ecosystem-modelling)
-  Instructions to use EcoSysEM platform via Command Line Interface (CLI) | [GO](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli)
-  Function Navigation | [GO](#function-navigation)
-  Error List | [GO](#error-list)
-  Contact | [GO](#contact)
____________________________

## Before having fun...
> [!NOTE]
> To open the links in a new tab: right click on the link + "Open link in new tab".

### :gear: Anaconda Python installation
EcoSysEM platform is built up in Python. To execute this Python scripts is recommended the installation of **Anaconda**. **Anaconda Python** is a free, open-source platform that allows to write and execute code in the programming language Python ([Python Tutorial](https://docs.python.org/3/tutorial/index.html)). This platform simplifies package installation, managment and development, and alos comes with a large number of libraries/packages that can be you for your projects. To install **Anaconda**, just head to the [Anaconda Documentation website](https://docs.anaconda.com/free/anaconda/install/index.html) and follow the instructions to download teh installer for your operating system.

[🔼 Back to **Contents**](#readme-contents)

### Anaconda Navigator
Anaconda Navigator is a desktop graphical user interface that allows you to launch applications and efficiently manage conda packages, environments, and channels without using command-line commands. For more info, click [here](https://docs.anaconda.com/free/navigator/).

[🔼 Back to **Contents**](#readme-contents)

### Anaconda Prompt or Terminal
Anaconda Prompt is a command line interface with Anaconda Distribution. Terminal is a command line interface that comes with macOS and Linux. To open it in **Windows**: Click Start, search for _"Anaconda Prompt"_ and click to open. In **macOS**: use Cmd+Space to open Spotlight Search and type _"Navigator"_ to open the program. In **Linux-CentOS**: open Applications > System Tools > Terminal.

[🔼 Back to **Contents**](#readme-contents)

### Spyder
Spyder is a Python development environment with many features for working with Python code, such as a text editor, debugger, profiler, and interactive console. You can launch **Spyder** using the **Anaconda Navigator**. For Spyder Tutorials, click [here](https://www.youtube.com/watch?v=E2Dap5SfXkI&list=PLPonohdiDqg9epClEcXoAPUiK0pN5eRoc&ab_channel=SpyderIDE).

[🔼 Back to **Contents**](#readme-contents)

### Python packages
A **Python package** is a collection of files containing Python code (i.e., modules). To execute **EcoSysEM platform**, the following packages must to be installed:
- **<ins>NumPy</ins>**. NumPy is the fundamental package for scientific computing in Python. It is a Python library that provides a multidimensional array object, various derived objects (such as masked arrays and matrices), and an assortment of routines for fast operations on arrays, including mathematical, logical, shape manipulation, sorting, selecting, I/O, discrete Fourier transforms, basic linear algebra, basic statistical operations, random simulation and much more. For more info and tutorials, click [here](https://numpy.org/).
- **<ins>Pandas</ins>**. Pandas is a fast, powerful, flexible and easy to use open source data analysis and manipulation tool, built on top of the Python programming language. For more info and tutorials, click [here](https://pandas.pydata.org/).
- **<ins>Matplotlib</ins>**. Matplotlib is a library for creatinc static, animated and interactive visualizations in Python. For more info and tutorials, click [here](https://matplotlib.org/).
- **<ins>SciPy</ins>**. SciPy is a collection of mathematical algorithms and convenience functions built on NumPy . It adds significant power to Python by providing the user with high-level commands and classes for manipulating and visualizing data. For more info and tutorials, click [here](https://scipy.github.io/devdocs/tutorial/index.html).
- **<ins>iteration_utilites</ins>**. Iteration_utilities is a collection of functional programming based on and utilizing iteratiors and generators. Most of the functions are inspiered by the _itertools_ module, but implemented in C to achieve a better overall performance. For more info and tutorials, click [here](https://iteration-utilities.readthedocs.io/).
- **<ins>Xarray</ins>**. Xarray introduces labels in the form of dimensions, coordinates and attributes on top of raw NumPy-like multidimensional arrays, which allows for a more intuitive, more concise, and less error-prone developer experience. xarray is better suited for more complex tasks that involve labeled arrays or multi-dimensional arrays with missing or incomplete data. For more info ant tutorials, click [here](https://docs.xarray.dev/en/stable/user-guide/index.html).
- **<ins>Earthaccess</ins>**. Earthaccess is a python library to search for, and download or stream NASA Earth science data with just a few lines of code. For more info and tutorials, click [here](https://earthaccess.readthedocs.io/en/latest/).
- **<ins>cdsapi</ins>**. The Climate Data Store (CDS) Application Program Interface (API) is a service providing programmatic access to CDS and ADS data. For more info and tutorials, click [here](https://cds.climate.copernicus.eu/how-to-api).
- **<ins>molmass</ins>**. Molmass is a Python library, console script, and web application to calculate the molecular mass (average, nominal, and isotopic pure), the elemental composition, and the mass distribution spectrum of a molecule given by its chemical formula, relative element weights, or sequence. For more info and tutorials, click [here](https://pypi.org/project/molmass/).
- **<ins>pyatmos</ins>**. Pyatmos is an archive of scientific routines that estimates the vertical structure of atmosphere with various atmospheric density models. For more info and tutorials, click [here](https://pypi.org/project/pyatmos/).

[🔼 Back to **Contents**](#readme-contents)

### Installation of packages using Anaconda Navigator
You can install any Python package using the **Anaconda Navigator**. For this, execute the navigator and click to **Environments**. In this section you can install new packages and delete the already installed. For more info, click [here](https://docs.anaconda.com/free/navigator/).

[🔼 Back to **Contents**](#readme-contents)

### Installation of packages using pip
**pip** is the package installer for Python. In general, pip installs the minimal instalation requirements automatically, but not the optionals requirements. To install the mentioned packages using pip, you have only to write the following command lines in **Anaconda Prompt or Terminal**:
#### · <ins>Anaconda Prompt</ins>
**NumPy**:
```
pip install numpy
```
**Pandas**:
```
pip install pandas
```
**Matplotlib**:
```
pip install matplotlib
```
**SciPy**:
```
pip install scipy
```
**iteration_utilites**:
```
pip install iteration_utilities
```
**Xarray**:
```
pip install xarray
```
**earthaccess**:
```
pip install earthaccess
```
**cdsapi**:
```
pip install cdsapi
```
**molmass**:
```
pip install -U "molmass[all]"
```
**pyatmos**:
```
pip install pyatmos
```
**netCDF4 (h5netcdf)**
```
pip install netCDF4 h5netcdf
```
#### · <ins>Windows Terminal</ins>
**NumPy**:
```
python -m pip install numpy
```
**Pandas**:
```
python -m pip install pandas
```
**Matplotlib**:
```
python -m pip install matplotlib
```
**SciPy**:
```
python -m pip install scipy
```
**iteration_utilites**:
```
python -m pip install iteration_utilities
```
**Xarray**:
```
python -m pip install xarray
```
**earthaccess**:
```
python -m pip install earthaccess
```
**cdsapi:**
```
python -m pip install cdsapi
```
**molmass**:
```
python -m pip install -U "molmass[all]"
```
**pyatmos**:
```
python -m pip install pyatmos
```
**netCDF4 (h5netcdf)**
```
python -m pip install netCDF4 h5netcdf
```
[🔼 Back to **Contents**](#readme-contents)

____________________________

## :clipboard: Instructions for downloading and setting up EcoSysEM platform
1. Download .zip code. Last version: `v0.2` **(Pre-release)**. [Download release](https://github.com/soundslikealloy/EcoSysEM/archive/refs/tags/v0.2.zip).
2. Extract files to a destination (Recommendation - Desktop).
3. Modify (if necessary) parameter databases using Excel files in folder  `\ecosysem\db\Excels` (see [How to modify parameter databases](#how-to-modify-parameter-databases) section).
4. Modify existing reaction databases or create a new one using Excel files in folder `\ecosysem\reactions\Excels` (see [How to modify reaction databases](#how-to-modify-reaction-databases) section).
5. Execute **EcoSysEM platform** via Spyder (see [Instructions to use EcoSysEM platform via Spyder](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) section) or Command Line Interface (see [Instructions to use EcoSysEM platform via Command Line Interface (CLI)](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli) section).

[🔼 Back to **Contents**](#readme-contents)

### How to modify parameter databases
All important parameters are saved using local databases in .csv format (read by code) and .xlsx format (to create/modify databases). For now, the following parameter databases are included:
- Standard Gibbs free energy of formation (ΔG<sub>f</sub><sup>0</sup>).
- Standard enthalpy of formation (ΔH<sub>f</sub><sup>0</sup>).
- Henry's law constant of solubility at standard temperature (H<sub>S</sub><sup>0</sup>).
- Temperature dependence of Henry's law constant of solubility (B).

To modify existing databases (_i.e.,_ include new parameter values), open de .xlsx file in `ecosysem\db\Excels\*.xlsx`. Once finished, **1)** _Save_ Excel file in `\ecosysem\db\Excels\` folder and **2)** _Save as_ the document in .csv format in `\ecosysem\db\` folder. All parameter databases have the same structure:

|  IUPAC  |  Formula  |  Phase  |  Value  |  REF  |
| ------- | --------- | ------- | ------- | ----- |
| IUPAC 1 | Formula 1 | Phase 1 | Value 1 | REF 1 |
| IUPAC 2 | Formula 2 | Phase 2 | Value 2 | REF 2 |

Where **IUPAC** is the name of compound using IUPAC nomenclautre, **Formula** is the chemical formula of compound (see [Formulization of compounds](#formulization-of-compounds)), **Phase** is the fluid phase in which parameter has been measured or estimated, **Value** is the value of parameter, and **REF** is the literature reference (see [database references](#database-references)).
The ΔG<sub>f</sub><sup>0</sup> and ΔH<sub>f</sub><sup>0</sup> parameters have three possible phases: G - gas, L - liquid, and S - solid. The H<sub>S</sub><sup>0</sup> and B parameters have two possible phases: FW - freshwater (liquid), and SW - sewater (liquid).

#### <ins>Database references</ins>
- Standard Gibbs free energy of formation (ΔG<sub>f</sub><sup>0</sup>).
    - Thaurer1997: R. Thauer, K. Jungermann & K. Decker (1977), doi: 10.1128/br.41.1.100-180.1977.
    - Kleerebezem2010: R. Kleerebezem & M. van Loosdrecht (2010), doi: 10.1080/10643380802000974.
    - Zumdahl2012: S. Zumdahl & D. DeCoste (2017), CHemical Principles (7th edition), ISBM: 9781305856745.
    - Beber2020: M. Beber _et al._ (2021), eQuilibrator 3.0, doi: 10.1093/nar/gkab1106.
- Standard enthalpy of formation (ΔH<sub>f</sub><sup>0</sup>).
    - Alberty2004: R. Alberty (2004), doi: 10.1016/j.bpc.2004.05.003.
    - Speight2005: J. Speight (2005), Lange's Handbook of Chemistry (16th edition), ISBM: 0-07-143220-5.
    - Kleerebezem2010: R. Kleerebezem & M. van Loosdrecht (2010), doi: 10.1080/10643380802000974.
- Henry's law constant of solubility at standard temperature (H<sub>S</sub><sup>0</sup>).
    - Murray1969: C. Murray, J. Riley & T. Wilson (1969), doi: 10.1016/0011-7471(69)90020-5.
    - Murray1970: C. Murray & J. Riley (1970), doi: 10.1016/0011-7471(70)90100-2.
    - Crozier1974: T. Crozier & S. Yamamoto (1974), doi: 10.1021/je60062a007.
    - Douabul1979: A. Douabul & J. Riley (1979), doi: 10.1021/je60083a014.
    - Douabul1979b: A. Douabul & J. Riley (1979), doi: 10.1016/0198-0149(79)90023-2.
    - Hedengren2000: D. Hedengren (2000), No. HNF-5174-FP.
    - Zacharia2005: I. Zacharia & W. Deen (2005), doi: 10.1007/s10439-005-8980-9.
    - Sander2023: R. Sander (2023), doi: 10.5194/acp-23-10901-2023.
- Temperature dependence of Henry's law constant of solubility (B).
    - Asm: Assumed equal to the value in freshwater.
    - Sander2023: R. Sander (2023), doi: 10.5194/acp-23-10901-2023.

[🔼 Back to **Instructions (Downloading & Setting up)**](#clipboard-instructions-for-downloading-and-setting-up-ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### How to modify reaction databases
Like parameter databases, all chemical and biotic reactions are defined using local databases in .csv format (read by code) and .xlsx format (to create/modify databases). In reaction databases, stocihiometric matrix of reactions is defined. For now, the following reaction databases are included:
- pHSpeciation. Stoichiometric matrix of acid-base equilibrium reactions.
- Metabolisms. Stoichiometrix matrix of metabolic reactions.

To modify existing databases open de .xlsx file in `ecosysem\reactions\Excels\*.xlsx`. The user can also create a new database reaction. Once finished, **1)** _Save_ Excel file in `\ecosysem\reactions\Excels\` folder and **2)** _Save as_ the document in .csv format in `\ecosysem\reactions\` folder. All reaction databases have the same structure:

|  Compounds  |  Reaction n<sup>o</sup> 1  |  Reaction n<sup>o</sup> 2  |  ...  |  Reaction n<sup>o</sup> N  |
| ----------- | -------------------------- | -------------------------- | ----- | -------------------------- |
| Compound 1  | Stoich. value 1.1 | Stoich. value 1.2 |  ...  | Stoich. value 1.N |
| Compound 2  | Stoich. value 2.1 | Stoich. value 2.2 |  ...  | Stoich. value 2.N |
| Compound 3  | Stoich. value 3.1 | Stoich. value 3.2 |  ...  | Stoich. value 3.N |

Where **Stoich. value A.B** is the stoichiometric value of _Compound B_ for _Reaction A_. Stoichiometric values are negative for substrates (<0) and positive for products (>0). If a compound does not participate in a reaction, that cell is left blank. Each column is a specific reaction, and in the headers is written the reaction name and its abbreviation in parentheses.  

· For example, biotic ammonia oxidation to nitrate by comammox bacteria (CMX): NH<sub>3</sub> + 2.0·O<sub>2</sub> → NO<sub>3</sub><sup>-</sup> + H<sub>2</sub>O + H<sup>+</sup>.

|  Compound  |  Complete ammonia oxidation (CMX)  |
| ---------- | ---------------------------------- |
| H+ | 1 |
| H2O | 1 |
| O2 | -2 |
| NH3 | -1 |
| NO2- |  |
| NO3- | 1 |

> [!NOTE]
> In pHSpeciation database, a different header format is used. The participating compounds are written between '/' and the type of protonation between parenthesis:
> - First deP: first deprotonation.
> - Second deP: second deprotonation.
> - Third deP: third deprotonation.
>
> · For example, for the acid-base equilibrium of sulfuric acid (H<sub>2</sub>SO<sub>4</sub>): H<sub>2</sub>SO<sub>4</sub> ⇔ HSO<sub>4</sub><sup>-</sup> + H<sup>+</sup> ⇔ SO<sub>4</sub><sup>2-</sup> + 2.0·H<sup>+</sup>.
> |  Compound  |  H2SO4/HSO4-/SO4-2 (First deP)  |  H2SO4/HSO4-/SO4-2 (Second deP)  |
> | ---------- | ------------------------------- | -------------------------------- |
> | H+ | 1 | 1 |
> | H2SO4 | -1 |  |
> | HSO4- | 1 | -1 |
> | SO4-2 |  | 1 |

[🔼 Back to **Instructions (Downloading & Setting up)**](#clipboard-instructions-for-downloading-and-setting-up-ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### Formulization of compounds
To identify the differents compounds across the **EcoSysEM platform** and associated databases, the molecular formula of compounds has been chosen. It is crucial that each compound used in the platform has its unique formulation. Make sure that the compound formula match in all databases and scripts.

· Guidelines for formulating compounds:
- Use the same chemical symbol of elements. For example: _H_ for hydrogen, _O_ for oxygen, _Ar_ for argon, _Na_ for sodium...
- If the compound is charged (like ion ammonium, NH<sub>4</sub><sup>+</sup>; or sulphate, SO<sub>4</sub><sup>2-</sup>), the charge sign (+ or -) is always precedes the number of charges of the ion, that is, in reverse of the normal formulizaton of ions. **Use only '+' and '-' to represent the charge of the ion.** For example: _SO4-2_ for sulphate (SO<sub>4</sub><sup>2-</sup>), _SO3-2_ for sulfite (SO<sub>3</sub><sup>2-</sup>) or _CO3-2_ for carbonate (CO<sub>3</sub><sup>2-</sup>.
- For free radicals, the compound formulization starts with 'rad' and followed by the symbol of element/compound. For example: _radOH_ for hydroxyl radical (OH·), _radBr_ for bromine radical (Br·), or _radCH3_ for methyl radical (CH<sub>3</sub>·).
- For organic compounds, you can use the chemical formula or an abbreviation. For example: _CH3COOH_ or _AcOH_ for acetic acid, _C3H7COOH_ or _ButyOH_ for butyric acid or _C6H1206_ or _Glc_ for glucose.

[🔼 Back to **Instructions (Downloading & Setting up)**](#clipboard-instructions-for-downloading-and-setting-up-ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

 ## :clipboard: Instructions to use EcoSysEM platform via Spyder
1. Launch **Spyder**. For Spyder Tutorials, click [here](https://www.youtube.com/watch?v=E2Dap5SfXkI&list=PLPonohdiDqg9epClEcXoAPUiK0pN5eRoc&ab_channel=SpyderIDE).
2. Set (at least) the following panes in Spyder (most are selected by default): `Files`, `Editor`, `IPython Console`, `Plots`, `Help`, `Historial`.
   From Spyder taskbar: <ins>V</ins>iew / Panes ▸.
3. Go to the **Code folder<sup>2</sup>** using the `Files` pane and open `ecosysem_spyder.py` file.
    &#09;<br><sup><sup>2</sup>Code folder: folder with `ecosysem_spyder.py` file (Folder: `EcoSysEM\ecosysem`). </sup>
5. Program your script. For user guide, click [here](#fundamentals-and-usage).
6. Run `ecosysem_spyder.py` script with Ctrl + Intro, F5 or Play symbol of _Run toolbar_.  

[🔼 Back to **Contents**](#readme-contents)

## EcoSysEM user guide
This section is an overview and explanation of important features of **EcoSysEM platform**.

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### EcoSysEM package layout
Important modules and how to import functions or classes from them are listed below. Classes names start with a capital letter, functions with a lower letter, and attributes with a dot (.) and lower letter:
```python
from ecosysem.module import function
from ecosysem.module import Class

ecosysem
  ├── ecosysem_cmd.py
  ├── ecosysem_spyder.py
  ├── ecosystem.py
  │      └── EcoA
  │           ├── ecoysDirEdges
  │           └── ecoysDiHypergraph
  ├── environments.py
  │      ├── Environment (this behaviours can also be called by model classes. E.g., ISA.getDGr(...))
  │      │    ├── loadData
  │      │    ├── getAttributeNames
  │      │    ├── combData
  │      │    └── getDGr
  │      ├── ISA
  │      │    ├── .altitude
  │      │    ├── .temperature
  │      │    ├── .pressure
  │      │    ├── .Pi
  │      │    ├── .Ci_G
  │      │    ├── .Ci_LFW
  │      │    ├── .Ci_LSW
  │      │    ├── .compounds
  │      │    ├── .compositions
  │      │    ├── .pH
  │      │    ├── .H2O
  │      │    ├── .layers
  │      │    ├── .resolution
  │      │    ├── .DGr (if ISA.getDGr() is called)
  │      │    ├── setComposition
  │      │    ├── plotTandP
  │      │    └── plotCompsProfilesISA
  │      ├── MERRA2 # Attributes can be different
  │      │    ├── .altitude
  │      │    ├── .temperature
  │      │    ├── .pressure
  │      │    ├── .H
  │      │    ├── .LR
  │      │    ├── .lat
  │      │    ├── .lon
  │      │    ├── .PS
  │      │    ├── .T2M
  │      │    ├── .TROPH
  │      │    ├── .TROPPB
  │      │    ├── .TROPT
  │      │    └── getDataMERRA2
  │      ├── ISAMERRA2 # Attributes can be different
  │      │    ├── .altitude
  │      │    ├── .temperature
  │      │    ├── .pressure
  │      │    ├── .Pi
  │      │    ├── .Ci_G
  │      │    ├── .Ci_LFW
  │      │    ├── .Ci_LSW
  │      │    ├── .lat
  │      │    ├── .lon
  │      │    ├── .compounds
  │      │    ├── .compositions
  │      │    ├── ISAMERRA2.DGr (if .getDGr() is called)
  │      │    └── # Dynamic attributes from MERRA2 (if keysAsAttributes = True)
  │      ├── CAMS # Attributes can be different
  │      │    ├── .altitude
  │      │    ├── .lat
  │      │    ├── .lon
  │      │    ├── .P_level
  │      │    ├── .CH4
  │      │    ├── .CO
  │      │    ├── .CO2
  │      │    └─ getDataCAMS
  │      │── CAMSMERRA2
  │      │    ├── .altitude
  │      │    ├── .temperature
  │      │    ├── .pressure
  │      │    ├── .Pi
  │      │    ├── .Ci_G
  │      │    ├── .MolPct_G
  │      │    ├── .Ci_LFW
  │      │    ├── .Ci_LSW
  │      │    ├── .lat
  │      │    ├── .lon
  │      │    ├── .compounds
  │      │    ├── CAMSMERRA2.DGr (if .getDGr() is called)
  │      │    └── # Dynamic attributes from MERRA2 (if keysAsAttributes = True)
  │      └── GWB
  │           ├── .temperature
  │           ├── .pH
  │           ├── .salinity
  │           ├── .Ci_L
  │           ├── .methods
  │           ├── .fluidType
  │           ├── .compounds
  │           └── GWB.DGr (if .getDGr() is called)
  ├── reactions.py
  │      ├── KinP
  │      │    └── getKinP
  │      ├── KinRates
  │      │    └── getRs
  │      └── Reactions
  │           ├── getRxn
  │           ├── getRxnByComp
  │           └── getRxnByName
  └── thermodynamics.py 
         ├── ThP
         │    ├── getThP
         │    ├── getDeltaG0r
         │    ├── getDeltaH0r
         │    ├── ionicStrength
         │    ├── activity
         │    └── getKeq
         ├── ThEq
         │    ├── solubilityHenry
         │    ├── pHSpeciation
         │    └── plotpHSpeciation
         └── ThSA
              ├── exportDeltaGr
              └── getDeltaGr
```

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### Fundamentals and usage
This section clarifies concepts, design decisions and technical details of this package. **EcoSystem platform** is constituted by four main units:
- Environment definition and instance calling | [GO](#environment-definition-and-instance-calling)
  - General functions (for all environmental models) | [GO](#Environment)
  - Atmosphere
    - Ideal Earth's atmosphere (International Standard Atmosphere, ISA) | [GO](#ISA)
    - Modern-Era Retrospective analysis for Research and Applications, Version 2 (MERRA-2) | [GO](#MERRA2)
    - ISA-MERRA2 atmospheric model | [GO](#ISAMERRA2)
    - Copernicus Atmosphere Monitoring Service (CAMS) | [GO](#CAMS)
    - CAMS-MERRA2 atmospheric model | [GO](#CAMSMERRA2)
  - Hydrosphere
    - General (or non-specific) Water Body (GWB) | [GO](#GWB)
  - How to create a new environment (class or subclass) | [GO](#create-new-environment)
- Thermodynamic State Analysis (ThSA) | [GO](#thermodynamic-state-analysis-thsa)
- Ecosystem modelling 🚧 | [GO](#ecosystem-modelling)
<!-- - Bio-Thermodynamic State Analysis (BioThSA) 🚧 | [GO](#bio-thermodynamic-state-analysis-biothsa)
- Ecosystem Analysis (EcoA) 🚧 | [GO](#ecosystem-analysis-ecoa) -->
[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Environment definition and instance calling</ins>
One of the advantages of Python is that supports both **Object-Oriented Programming (OOP)** and functional programming paradigms. The definition of environments are based on <ins>OOP paradigm</ins>. OOP is based on the following four principles: _Encapsulation_, _Inheritance_, _Abstraction_ and _Polymorphism_. For more info about OOP principles, click [here](https://fluxtech.me/blog/object-oriented-programming-vs-functional-programming/). 
- _Encapsulation_ principle allows to hide he internal state and behaviour of an object, and the object can only be acessed through a well-defined interface. With this, the user can change the properties of the environment without affecting the code hat uses the object.
- _Inheritance_ principle allows a new class to be defined based on an existing class, inheriting its attributes and methods.
- _Abstraction_ principle makes possible to work with objecs of a class without knowing the details of their implementation, which can make the code more robust and less error-prone.
- _Polymorphism_ principle enables the use of a common interface for different classes, making it possible to write code that can work with object of different types without knowing their specific class.

The benefits of OOP are _i_) organization, _ii_) state definition and tracking, _iii_) encapsulation of proceudre and data (_i.e.,_ specific functions and data can be stored together in a single class), _iv_) inheritance (making development more efficient and easier to maintain). For more information about OOP in Python, click [here](https://realpython.com/python3-object-oriented-programming/).

#

<a name="Environment">**General functions for all environmental models**</a><br>
Each environment has their own models in object form, with the corresponding attributes and behaviours (i.e., functions). Some behaviours are shared between the distinct environmental classes, which they are gathered in _Environment_ object. _Environment_ object (parent class) has a set of inherited classes (child class), where the latters inherits all attributes and behaviours of parent class. The child classes of _Environment_ are the distinct regions of Earth: _Atmosphere_, _Hydrosphere_, _Cryosphere_ or _Lithosphere_. For now, only _Atmosphere_ is available. At the same time, region objects (like _Atmosphere_) have their own inheritance, the distinct environmental models. For example, _Atmosphere_ object has _ISA_, _MERRA2_, _CAMS_, _ISAMERRA2_ and _CAMSMERRA2_. Inheritance is represented as `Child(Parent)`:
```python
Environment
  ├── Atmosphere(Environment)
  │      ├── ISA(Atmosphere)
  │      ├── MERRA2(Atmosphere)
  │      ├── CAMS(Atmosphere)
  │      ├── ISAMERRA2(Atmosphere)
  │      └── CAMSMERRA2(Atmosphere)
  ├── Hydrosphere(Environment)
  │      └── GWB(Hydrosphere)
  ├── Cryosphere(Environment) # As example. Currently unavailable.
  │      └── Glacier(Cryosphere)
  └── Lithosphere(Environment) # As example. Currently unavailable.
         └── Crust(Lithosphere)
```

### Environment.loadData &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.loadData(model, dataType, y, m=None, d=None, keys='All')
```
Get data in dictionary form.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **model : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Environmental model.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ ('dly', 'mly', 'cmly', 'yly' 'cyly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'dly' - Daly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Combined monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'yly' - Annual data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cyly' - Combined annual data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ ('mly' and 'dly') or _list of int_ ('cmly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_, _optional, default; None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictVar : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dictionary with requested variables.<br>

### Environment.combData &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.combData(model, dataType, year, month, days=None, keys='All', dataDelete=False)
```
Get average and standard deviation from a group of data.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **model : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Environmental model.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ ('mly', 'cmly', 'yly', 'cyly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Generate monthly data from daily data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Generate combined monthly data from monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'yly' - Generate annual data from monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cyly' - Generate combined annual data from monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **year : _int_ ('mly' and 'dly') or _list of int_ ('cmly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **month : _int_ ('mly' and 'dly') or _list of int_ ('cmly', 'yly', 'cyly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **days : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataDelete : _bool_, _optional, default: False_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Delete daily or monthly data after the average calculation.<p> 
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **NPZ file in folders, `data\{selected environment}\{selected model}\mly\`, `data\{selected environment}\{selected model}\cmly\`, `data\{selected environment}\{selected model}\yly\` or `data\{selected environment}\{selected model}\cyly\`**<br>

### Environment.getAttributeNames &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.getAttributeNmes()
```
Return attribute names of an Environment object as a list.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **None** <br>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **attributes : _list of str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Attribute names of Environment object.<br>

### Environment.getDGr &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.getDGr(typeRxn, input_, specComp)
```
Compute (non-)stadard Gibbs free energy using the information from environmental models (e.g., temperature, pH, concentrations, and so on).<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; What reaction database is used, matching with CSV name in `reactions\` folder.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **input_ : _str_ or _list of strs_ or _ndarray of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested compound(s) or reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total concentrations of compounds `{'compounds': [concentrations]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **specComp : _bool_, _str_, _list of strs_, _ndarray of strs_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Bool_ if `input_` are compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _str_, _list of strs_ or _ndarray of strs_ if `input_` are reactions.<p>
**Returns:** <br>
**New attribute (`.DGr`) is created in object instance.**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.DGr : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gibbs free energy values. `{'rxnName_pH:#.#': [DGr]}`<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

<a name="ISA">**Ideal Earth's atmosphere (International Standard Atmosphere, ISA)**</a><br>
The International Standard Atmosphere (ISA) is a static atmospheric model of how pressure, temperature, density and viscosity of the Earth's atmosphere change over a a wide range of altitudes. The ISA model is detailed in ISO 2533:1975[^1]. The ISA model divides the atmosphere into different layers with specific physical properties (such as temperature rate change, base temperature, base atmospheric pressure or base atmospheric density). With these properties, the temperature and pressure profiles are calculated. On the other hand, the ISA model assumes that Earth's atmosphere is a mixture of gas, water vapor and a certain quantity of aerosols, in which the composition remains pratically constant up to altitudes of 90 - 95 km. The dry (0%<sub>vol</sub> of water) and wet composition (up to 4%<sub>vol</sub> of water) of Earth's atmosphere are from National Oceanic and Atmospheric Administration (NOAA)[^2], Agency for Toxic Substances and Disease Registry (ATSDR), United States ENvironmental Protection Agency, _Atmospheric Radiation: Theoretical Basis (2<sup>nd</sup> edition)_[^3], and scientific publications[^4]  (**Table 1**).

**Table 1. Dry Earth's atmosphere composition (in %<sub>vol</sub>)**

| Compound | Value | Compound | Value | Compound | Value |
|---|---|---|---|---|---|
| __N<sub>2</sub>__  | 7.8084·10<sup>-1</sup> | __H<sub>2</sub>__  | 5.500·10<sup>-7</sup> | __NH<sub>3</sub>__  | 6.000·10<sup>-9</sup>  |
| __O<sub>2</sub>__  | 2.0946·10<sup>-1</sup> | __N<sub>2</sub>O__ | 3.300·10<sup>-7</sup> | __HNO<sub>2</sub>__ | 1.000·10<sup>-9</sup>  |
| __Ar__             | 9.3400·10<sup>-3</sup> | __CO__             | 1.000·10<sup>-7</sup> | __HNO<sub>3</sub>__ | 1.000·10<sup>-9</sup>  |
| __CO<sub>2</sub>__ | 4.2600·10<sup>-4</sup> | __Xe__             | 9.000·10<sup>-8</sup> | __H<sub>2</sub>S__  | 3.300·10<sup>-10</sup> |
| __Ne__             | 1.8182·10<sup>-5</sup> | __O<sub>3</sub>__  | 7.000·10<sup>-8</sup> | | |
| __He__             | 5.2400·10<sup>-6</sup> | __NO<sub>2</sub>__ | 2.000·10<sup>-8</sup> | | |
| __CH<sub>4</sub>__ | 1.9200·10<sup>-6</sup> | __SO<sub>2</sub>__ | 1.500·10<sup>-8</sup> | | |
| __Kr__             | 1.1400·10<sup>-6</sup> | __I<sub>2</sub>__  | 1.000·10<sup>-8</sup> | | |

### ISA &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
instance_ISA = ISA(layers='All', phase='All', H2O=0.0, pH=7.0, selCompounds=None, selAlt=None, resolution=1000, showMessage=True)
```
Create an instance of `ISA` object.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **layers : _str_ ('All'), _int_ or _list of int_ (from 0 to 7), _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Selection of atmosphere layers defined by ISA model[^1].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Selection of phase of vertical profile composition.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'G' - Gas phase.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L' - Liquid phase.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-FW' - Liquid freshwater. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-SW' - Liquid seawater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'All' - All phases.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **H2O : _float_ (0.0 to 0.04), _optional, default: 0.0_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Water content of atmosphere.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _optional, default: 7.0_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; pH of aerosols.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **selCompounds : _str_ or _list of str_, _optional, default: None_ (i.e., all compounds are considered)**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Interested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **selAlt : _float_ or _list of float_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Selected altitude region (in meters). 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; When selAlt : _float_, the region is from 0.0 to selAlt.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; When selAlt : _list of float_ (length = 2), the region is [min_selAlt, max_selAlt].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **resolution : _int_ or _float_, _optional, default: 1000_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Resolution of altitude array, that is, the size of altitude nodes per layer (in meters).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **showMessage : _bool_, _optional, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Boolean to set whether informative messages are displayed in Console.<p>
**Attributes:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.altitude : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric altitudes in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.temperature : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric temperatures in Kelvin.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.pressure : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric pressures in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Pi : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Partial pressure of compounds throughout the atmosphere. `{'compound': [partial pressure]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_G : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gas concentration of compounds throughout the atmosphere. `{'compound': [gas concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_LFW : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of compounds in freshwater aerosol throughout the atmosphere. `{'compound': [liquid concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_LSW : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of compounds in seawater aerosol throughout the atmosphere. `{'compound': [liquid concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.compounds : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.compositions : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric dry composition of the atmospehre in %vol. `{'compound': [air composition]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.pH : _float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; pH of aerosols.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.H2O : _float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric water content in %vol.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.H2O : _float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric water content in %vol.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.layers : _str_, _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Selected layers of the atmosphere.<br>

Here is an example:
```python
from envdef import ISA
import numpy as np

newISA = ISA(0, resolution = 500)
>>> print(newISA.altitude)
[    0.   500.  1000.  1500.  2000.  2500.  3000.  3500.  4000.  4500.
  5000.  5500.  6000.  6500.  7000.  7500.  8000.  8500.  9000.  9500.
 10000. 10500. 11000.]
>>> print(newISA.temperature)
[288.15 284.9 281.65 278.4 275.15 271.9 268.65 265.4 262.15 258.9
 255.65 252.4 249.15 245.9 242.65 239.4 236.15 232.9 229.65 226.4
 223.15 219.9 216.65]

newISA = ISA(0, resolution = 500, selAlt = [1000, 3500])
>>> print(newISA.altitude)
[1000. 1500. 2000. 2500. 3000. 3500.]
>>> print(newISA.temperature)
[281.65 278.4 275.15 271.9 268.65 265.4]

C = newISA.Ci_LSW
# np.round(a, decimals) -> NumPy function: Evenly round to the given number of decimals.
>>> print(np.round(C['O2'], 7)) 
[0.0003134 0.0003138 0.0003144 0.0003153 0.0003164 0.0003178]
```

### ISA.setComposition &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.setComposition(compound, composition)
```
Add (if _compound_ does not exist) or modify (if _compound_ does exist) composition of the `Environment` instance.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compound : _str_ or _list of strs_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; New compound or compound to modify.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **composition : _float_ or _list of float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Composition of new or existing compound.<br>

### ISA.plotTandP &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.plotTandP()
```
Plot temperature and pressure profiles of `ISA` instance.<p>
**Parameters:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **None** <p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Spyder plot** <br>

### ISA.plotCompsProfilesISA &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.plotCompsProfilesISA(Conc, xLabel, logCLabel=False, compounds=None)
```
Return vertical profiles in _format=dict_ of selected compounds or all compounds of `ISA` subclass.<p>
**Parameters:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **C : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Compound concenrations to be plotted.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **xLabel : _ndarray_ or _list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Altitude vector.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **logCLabel : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set logarithmic scale in x-coordinate. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compound : _str_ or _list of strs_, _optional, default: None_ (all compounds)** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The desired compound(s) to plot.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Spyder plot**<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

<a name="MERRA2">**Modern-Era Retrospective analysis for Research and Applications, Version 2 (MERRA-2)**</a><br>
The Modern-Era Retrospective analysis for Research and Applications, Version 2 (MERRA-2) provides data beginning in 1980. It was introduced to replace the original MERRA dataset because of the advances made in the assimilation system that enable assimilation of modern hyperspectral radiance and microwave observations, along with GPS-Radio Occultation datasets. It also uses NASA's ozone profile observations that began in late 2004. Additional advances in both the GEOS model and the GSI assimilation system are included in MERRA-2. Spatial resolution remains about the same (about 50 km in the latitudinal direction) as in MERRA.
Along with the enhancements in the meteorological assimilation, MERRA-2 takes some significant steps towards GMAO’s target of an Earth System reanalysis. MERRA-2 is the first long-term global reanalysis to assimilate space-based observations of aerosols and represent their interactions with other physical processes in the climate system. MERRA-2 includes a representation of ice sheets over (say) Greenland and Antarctica.

> [!NOTE]
> The user needs an **Earthdata** profile to access MERRA-2 data ([registration and login page](https://urs.earthdata.nasa.gov/home)). After registration, [login](https://urs.earthdata.nasa.gov/home) with your account and accept all _end-user license agreements_ (EULAs). You can find them in <ins>EULAs</ins> -> <ins>Accept New EULAs</ins>. When `MERRA2.getDataMERRA2()` is executed, the Earthdata login is requested with  the following `input()` lines:
> ```phytion
> Enter your Earthdata Login username:
> Enter your Earthdata password:
> ```

### MERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
instance_MERRA2 = MERRA2(dataType=None, y=None, m=None, d=None, bbox=(-180, -90, 180, 90), keys='All', keysAsAttributes=True, altArray=None, numAlt=50, showMessage=True)
```
Create an instance of `MERRA2` object.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'dly' - Daily data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'yly' - Yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Combined monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cyly' - Combined yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ or _list of int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _int_ for dataType 'dly', 'mly' and 'yly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _list of int_ [start_year, end_year] for dataType 'cmly' and 'cyly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _tuple_, _optional_, default: (-180, -90, 180, 90)_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earth's region of data, the bounding box.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _str_ or _list of str_, _optional_, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keysAsAttributes : _bool_, _optional_, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set if keys from MERRA2 database are saved as object attributes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **altArray : _list of float_ or _np.ndarray of floats_, _optional_, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested altitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **numAlt : int_, _optional_, default: 50_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of altitude steps from 0.0 m to maximum tropopause altitude, if list of altitude is not given by the user with `altArray`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **showMessage : _bool_, _optional, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Boolean to set whether informative messages are displayed in Console.<p>
**Attributes:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.altitude : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric altitudes in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.temperature : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric temperatures in Kelvin.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.pressure : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric pressures in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lat : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of latitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lon : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of longitudes.<p>
**Default dynamic attributes**:<br> 
It can be other variables from MERRA2 databases (see [MERRA-2 documentation](https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf)).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.H : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earth's topography in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.LR : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Lapse rate distribution in K/km.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.PS : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Surface pressure in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.T2M : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2-meter air temperature in Kelvins, considered as surface temperature.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.TROPH : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tropopause altitude in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.TROPPB : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tropopause pressure in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.TROPT : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tropopause temperature in Kelvins.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.LR_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of lapse rate in K/km.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.PS_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of surface pressure in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.T2M_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of surface temperature in Kelvins.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.TROPH_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of tropopause altitude in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.TROPPB_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of tropopause pressure in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.TROPT_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of tropopause temperature in Kelvins.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<p>

The _MERRA2_ object has two modes: **Downloading** and **Loading/Combining**. To create a new _MERRA2_ instance in **Downloading** mode, this must be created without arguments. The data from MERRA-2 can be downloaded using `MERRA2.getDataMERRA2`. The data will be saved in the folder `data\MERRA2` in .npz format (more info [here](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html)). Here is an example:
```python
from envdef import MERRA2

newMERRA2 = MERRA2()

# Get monthly data from online databases
## (Default arguments: product = 'M2I1NXASM', version = '5.12.4', var = ['PS', 'T2M', 'TROPT', 'TROPPB'])
newMERRA2.getDataMERRA2(dataType = 'mly', years = 1995, months = 10, days = 'All')
```
Once the data is downladed, the user can load the data creating a new _MERRA2_ instance in **Loading** mode with, at least, `dataType` and `y` arguments. The user can get the data calling the attributes defined above. Here is an example:
```python
from envdef import MERRA2

# Yearly data from 2020 was previously downloaded
newMERRA2 = MERRA2(dataType = 'yly', y = 2020, bbox = (-180, -90, -178.125, -88.5))

>>> print(newMERRA2.getAttributeNames())
['environment', 'model', 'mode', 'dataType', 'bbox', 'temperature', 'pressure', 'altitude', 'H',
'PS', 'PS_std', 'TROPPB', 'TROPPB_std', 'T2M', 'T2M_std', 'TROPT', 'TROPT_std', 'TROPH', 'TROPH_std',
'LR', 'LR_std', 'lat', 'lon']

# Get data
>>> print(newMERRA2.lat)
[-90.  -89.5 -89.  -88.5]
>>> print(newMERRA2.lon)
[-180.    -179.375 -178.75  -178.125]
>>> print(newMERRA2.PS)
[[68299.766 68299.766 68299.766 68299.766]
 [67556.914 67563.96  67571.164 67578.46 ]
 [66478.73  66489.18  66499.94  66510.99 ]
 [65456.55  65458.25  65464.496 65470.754]]
>>> print(newMERRA2.temperature)
# newMERRA2.temperature has NaN values because altitude < surface height ASL (topography).
[[[         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]]
 [[         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]]
...
 [[208.89917863 208.89917863 208.89917863 208.89917863]
  [208.91566434 208.91529312 208.91491576 208.9145402 ]
  [209.0320008  209.03220517 209.0325211  209.03274575]
  [209.21578416 209.21842957 209.21895037 209.21945189]]
 [[208.4474782  208.4474782  208.4474782  208.4474782 ]
  [208.48575225 208.48551779 208.48527922 208.48503976]
  [208.59335294 208.59345043 208.59367159 208.59381527]
  [208.72229302 208.72450837 208.72443997 208.72438963]]]
```

### MERRA2.getDataMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>

> [!NOTE]
> If a significant amount of data needs to be downloaded, we recommend to run MERRA2.getDataMERRA2 [via Command Line Interface](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli).
> With this, you can download data sets in parallel - one set (_e.g._, one entire month) per Command Prompt.

```python
MERRA2.getDataMERRA2(dataType, years, months, days='All', product='M2I1NXASM', version='5.12.4', bbox=(-180, -90, 180, 90), var=['PS', 'T2M', 'TROPT', 'TROPPB'])
```
Download data from MERRA2 database.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ or _list of int_ ('dly', 'mly', 'cmly', 'All')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'dly' - Daly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Combined monthly data <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'All' - ['dly', 'mly', 'cmly'].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **years : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **months : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **days : _int_, _list of int_ or _str ('All')_, _optional, default: 'All'**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **product : _str_, _optional, default: 'M2I1NXASM'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Product of data (section of MERRA2 database).<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **version : _str_, _optional, default: '5.12.4'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Version of data.<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _tuple_, _optional, default: (-180, -90, 180, 90)_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earths region of data, the bounding box.<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **var : _list of str_, _optional, default: ['PS', 'T2M', 'TROPT', 'TROPPB']_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **NPZ file in folders, `data\MERRA2\dly\`, `data\MERRA2\mly\` and/or `data\MERRA\cmly\`**<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

<a name="ISAMERRA2">**ISA-MERRA2 atmospheric model**</a><br>
Combination of International Standard Atmosphere ([ISA](#ISA)) model and Modern-Era Retrospective analysis for Research and Applications Version 2 ([MERRA2](#MERRA2)).

### ISAMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
instance_ISAMERRA2 = ISAMERRA2(dataType, y, m=None, d=None, bbox=(-180, -90, 180, 90), compound=None, phase='All', altArray=None, numAlt=50, surftrop=None, keysAsAttributes=False, showMessage=True)
```
Create an instance of `ISAMERRA2` object.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'dly' - Daily data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'yly' - Yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Combined monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cyly' - Combined yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _int_ for dataType 'dly', 'mly' and 'yly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _list of int_ [start_year, end_year] for dataType 'cmly' and 'cyly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _tuple_, _optional_, default: (-180, -90, 180, 90)_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earth's region of data, the bounding box.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compound : _str_ or _list of str_, _optional_, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Interested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Selection of phase of vertical profile composition.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'G' - Gas phase.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L' - Liquid phase.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-FW' - Liquid freshwater. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-SW' - Liquid seawater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'All' - All phases.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **altArray : _list of float_ or _np.ndarray of floats_, _optional_, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested altitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **numAlt : int_, _optional_, default: 50_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of altitude steps from 0.0 m to maximum tropopause altitude, if list of altitude is not given by the user with `altArray`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **surftrop : _str_ ('surface', 'tropopause'), _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Get concentration from 2-meters air following topography (`surftrop='surface'`) or tropopause height (`surftrop='tropopause'`).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keysAsAttributes : _bool_, _optional_, default: False_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set if keys from MERRA2 database are saved as object attributes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **showMessage : _bool_, _optional, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Boolean to set whether informative messages are displayed in Console.<p>
**Attributes:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.altitude : _list_ or _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric altitudes in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If surftrop = None, `.altitude` is a list of altitudes in meters from 0.0 to max tropopause altitude.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If surftrop = 'surfase', `.altitude` is the Earth's topography in meters. Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If surftrop = 'tropopause', `.altitude` is the tropopause altitude in meters. Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.temperature : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric temperatures in Kelvin.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.pressure : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric pressures in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Pi : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Partial pressure of compounds throughout the atmosphere in Pascals. `{'compound': [partial pressure]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_G : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gas concentration of compounds throughout the atmosphere in mol/L. `{'compound': [gas concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_LFW : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of compounds in freshwater aerosol throughout the atmosphere in mol/L. `{'compound': [liquid concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_LSW : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of compounds in seawater aerosol throughout the atmosphere in mol/L. `{'compound': [liquid concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lat : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of latitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lon : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of longitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.compounds : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.compositions : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric dry composition of the atmospehre in %vol. `{'compound': [air composition]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<p>

If _MERRA2_ data has already been downloaded, the user can create a new _ISAMERRA2_ object with, at least, `dataType` and `y` arguments. The user can get the data calling the attributes defined above. Here is an example:
```python
from envdef import ISAMERRA2

# Yearly data from 2020 was previously downloaded
newISAMERRA2 = ISAMERRA2(dataType = 'yly', y = 2020, bbox = (-180, -90, -178.125, -88.5))
>>> print(newISAMERRA2.getAttributeNames())
['environment', 'compositions', 'compounds', 'model', 'temperature', 'pressure', 'altitude',
'Pi', 'Ci_G', 'Ci_LFW', 'Ci_LSW', 'lat', 'lon']

newISAMERRA2 = ISAMERRA2(dataType = 'yly', y = 2020, bbox = (-180, -90, -178.125, -88.5), keysAsAttributes = True)
>>> print(newISAMERRA2.getAttributeNames())
['environment', 'compositions', 'compounds', 'model', 'temperature', 'pressure', 'altitude',
'Pi', 'Ci_G', 'Ci_LFW', 'Ci_LSW', 'lat', 'lon', 'H', 'PS', 'PS_std', 'TROPPB', 'TROPPB_std',
'T2M', 'T2M_std', 'TROPT', 'TROPT_std', 'TROPH', 'TROPH_std', 'LR', 'LR_std', 'lat', 'lon']

# Get data
>>> print(newISAMERRA2.lat)
[-90.  -89.5 -89.  -88.5]
>>> print(newISAMERRA2.lon)
[-180.    -179.375 -178.75  -178.125]
>>> print(newISAMERRA2.T2M)
[[224.31621 224.31621 224.31621 224.31621]
 [223.43001 223.42714 223.4242  223.42134]
 [223.60979 223.6157  223.62134 223.62642]
 [225.37038 225.3871  225.40935 225.43036]]
>>> print(newISAMERRA2.temperature)
# newISAMERRA2.temperature has NaN values because altitude < surface height ASL (topography).
[[[         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]]
 [[         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]]
...
 [[208.89917863 208.89917863 208.89917863 208.89917863]
  [208.91566434 208.91529312 208.91491576 208.9145402 ]
  [209.0320008  209.03220517 209.0325211  209.03274575]
  [209.21578416 209.21842957 209.21895037 209.21945189]]
 [[208.4474782  208.4474782  208.4474782  208.4474782 ]
  [208.48575225 208.48551779 208.48527922 208.48503976]
  [208.59335294 208.59345043 208.59367159 208.59381527]
  [208.72229302 208.72450837 208.72443997 208.72438963]]]
>>> print(newISAMERRA2.Ci_LFW['O2'])
# newISAMERRA2.Ci_LFW['O2'] has NaN values because altitude < surface height ASL (topography).
[[[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
 [[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
...
 [[0.0005355  0.0005355  0.0005355  0.0005355 ]
  [0.00053454 0.00053452 0.00053449 0.00053447]
  [0.00053324 0.00053324 0.00053325 0.00053325]
  [0.00053258 0.00053257 0.00053256 0.00053254]]
 [[0.00052614 0.00052614 0.00052614 0.00052614]
  [0.00052481 0.00052478 0.00052475 0.00052472]
  [0.00052368 0.00052369 0.00052369 0.0005237 ]
  [0.00052402 0.00052402 0.00052402 0.00052401]]]
>>> print(newISAMERRA2.Ci_LFW['CH4'])
# newISAMERRA2.Ci_LFW['CH4'] has NaN values because altitude < surface height ASL (topography).
[[[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
 [[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
...
 [[6.10068828e-09 6.10068828e-09 6.10068828e-09 6.10068828e-09]
  [6.08954089e-09 6.08924764e-09 6.08896855e-09 6.08869882e-09]
  [6.07306869e-09 6.07310833e-09 6.07314832e-09 6.07322615e-09]
  [6.06301325e-09 6.06286162e-09 6.06269517e-09 6.06251621e-09]]
 [[6.00027347e-09 6.00027347e-09 6.00027347e-09 6.00027347e-09]
  [5.98451775e-09 5.98419948e-09 5.98389468e-09 5.98359967e-09]
  [5.97024831e-09 5.97031069e-09 5.97037078e-09 5.97046504e-09]
  [5.97235672e-09 5.97230126e-09 5.97226620e-09 5.97221065e-09]]]
```

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

<a name="CAMS">**Copernicus Atmosphere Monitoring Service (CAMS)**</a><br>
The Copernicus Atmosphere Monitoring Service (CAMS) provides continuous data and information on atmospheric composition, supporting a wide range of applications from air quality monitoring and forecasting to climate change assessment. It combines observations from satellites and in-situ measurements with sophisticated numerical models to provide consistent and quality-controlled data records. CAMS data typically encompasses global and regional analyses and forecasts of reactive gases, greenhouse gases (GHGs), aerosols, and stratospheric ozone. The data records span from approximately 2003 onwards for analyses and reanalyses, with forecast data available on a rolling basis.

> [!NOTE]
> The user needs an **Climate Data Store (CDS)** profile to access CAMS data ([registration and login page](https://cds.climate.copernicus.eu/)). After registration, [login](https://cds.climate.copernicus.eu/) with your account and accept Terms & conditions and Dataset licence, Licence to use Copernicus Products. You can find them in <ins>Your Profile</ins> -> <ins>Licences</ins>. The user also needs to add the installation and Script folder path in PATH using **set** (temporary) or **setx** (permanent) in a Command Prompt window (for **Windows**).
> Examples are available [here](https://docs.python.org/3/using/windows.html#setting-envvars).<br>
> Also, the user needs to copy a 2 line code and then paste into a  %USERPROFILE%\.cdsapirc file, where in the windows environment. <br>
> Please check for further help (also for **Linux** and **MacOS** users): [CDSAPI setup](https://ads.atmosphere.copernicus.eu/how-to-api).

### CAMS &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
instance_CAMS = CAMS(dataType=None, y=None, m=None, d=None, bbox=(-180, -90, 180, 90), keys='All', keys_to_reshape=['lat', 'lon', 'CH4', 'CH4_std', 'CO', 'CO_std', 'CO2', 'CO2_std'], keysAsAttributes=True, showMessage=True)
```
Create an instance of `MERRA2` object.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'dly' - Daily data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'yly' - Yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Combined monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cyly' - Combined yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _int_ for dataType 'dly', 'mly' and 'yly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _list of int_ [start_year, end_year] for dataType 'cmly' and 'cyly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _tuple_, _optional_, default: (-180, -90, 180, 90)_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earth's region of data, the bounding box.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_, _optional_, default: ['lat', 'lon', 'CH4', 'CH4_std', 'CO', 'CO_std', 'CO2', 'CO2_std']_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys_to_reshape : _str_ or _list of str_, _optional_, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of variables to be reshaped with ascending latitude and longitude order.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keysAsAttributes : _bool_, _optional_, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set if keys from CAMS database are saved as object attributes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **showMessage : _bool_, _optional, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Boolean to set whether informative messages are displayed in Console.<p>
**Attributes:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.altitude : _list_ or _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric altitudes in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lat : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of latitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lon : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of longitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.P_level : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of pressure levels in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.CH4 : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distribution of atmospheric methane (CH4) composition in kg/kg-air.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.CO : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distribution of atmospheric carbon monoxide (CO) composition in kg/kg-air.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.CO2 : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distribution of atmospheric carbon dioxide (CO2) composition in kg/kg-air.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.CH4_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of atmospheric methane (CH4) composition in kg/kg-air.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.CO_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of atmospheric carbon monoxide (CO) composition in kg/kg-air.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.CO2_std : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Standard deviation of atmospheric carbon dioxide (CO2) composition in kg/kg-air.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (latitude, longitude).<p>

Like _MERRA_ object, the _CAMS_ object has two modes: **Downloading** and **Loading/Combining**. To create a new _CAMS_ instance in **Downloading** mode, this must be created without arguments. The data from CAMS can be downloaded using `CAMS.getDataCAMS`. The data will be saved in the folder `data\CAMS` in .npz format (more info [here](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html)). Here is an example:
```python
from envdef import CAMS

# Get monthly data from online databases
## (Default arguments: pressure_levels = [50, 100, 200, 400, 600, 800, 900, 1000], variables = ["carbon_dioxide", "carbon_monoxide", "methane"])
newCAMS.getDataCAMS(dataType = 'mly', years = 2024, months = [4, 5], days = 'All', bbox = [90, -180, -90, 180])
```
Once the data is downladed, the user can load the data creating a new _CAMS_ instance in **Loading** mode with, at least, `dataType` and `y` arguments. The user can get the data calling the attributes defined above. Here is an example:
```python
from envdef import CAMS

# Yearly data from 2020 was previously downloaded
newCAMS = CAMS(dataType = 'yly', y = 2020, bbox = (-180, -90, -175, -85))
>>> print(newCAMS.getAttributeNames())
['base_path', 'environment', 'model', 'mode', 'dataType', 'bbox', 'lat', 'lon', 'P_level',
'altitude', 'CO2', 'CO2_std', 'CH4', 'CH4_std', 'CO', 'CO_std']

# Get data
>>> print(newCAMS.lat)
[-90. -89.6 -89.2 -88.8 -88.4 -88. -87.6 -87.2 -86.8 -86.4 -86. -85.6 -85.2 -84.8]
>>> print(newMERRA2.lon)
[-180. -179.6 -179.2 -178.8 -178.4 -178. -177.6 -177.2 -176.8 -176.4 -176. -175.6 -175.2]
>>> print(newISAMERRA2.temperature)
[[[9.91609274e-07 9.91609274e-07 9.91609274e-07 ... 9.91609274e-07
   9.91609274e-07 9.91609274e-07]
  [9.91629899e-07 9.91630003e-07 9.91630069e-07 ... 9.91629662e-07
   9.91629785e-07 9.91629832e-07]
  [9.91595073e-07 9.91595319e-07 9.91595509e-07 ... 9.91594571e-07
   9.91594760e-07 9.91594931e-07]
  ...
  [9.91185904e-07 9.91188652e-07 9.91191361e-07 ... 9.91177719e-07
   9.91180315e-07 9.91183100e-07]
  [9.91349755e-07 9.91352673e-07 9.91357666e-07 ... 9.91350551e-07
   9.91348173e-07 9.91346601e-07]
  [9.91458194e-07 9.91480221e-07 9.91502247e-07 ... 9.91448682e-07
   9.91451865e-07 9.91455029e-07]]
 [[9.91609274e-07 9.91609274e-07 9.91609274e-07 ... 9.91609274e-07
   9.91609274e-07 9.91609274e-07]
  [9.91629899e-07 9.91630003e-07 9.91630069e-07 ... 9.91629662e-07
   9.91629785e-07 9.91629832e-07]
  [9.91595073e-07 9.91595319e-07 9.91595509e-07 ... 9.91594571e-07
   9.91594760e-07 9.91594931e-07]
  ...
  [9.91185904e-07 9.91188652e-07 9.91191361e-07 ... 9.91177719e-07
   9.91180315e-07 9.91183100e-07]
  [9.91349755e-07 9.91352673e-07 9.91357666e-07 ... 9.91350551e-07
   9.91348173e-07 9.91346601e-07]
  [9.91458194e-07 9.91480221e-07 9.91502247e-07 ... 9.91448682e-07
   9.91451865e-07 9.91455029e-07]]
...
 [[7.61766302e-07 7.61766302e-07 7.61766302e-07 ... 7.61766302e-07
   7.61766302e-07 7.61766302e-07]
  [7.61136903e-07 7.61132360e-07 7.61127817e-07 ... 7.61150597e-07
   7.61146026e-07 7.61141450e-07]
  [7.60816723e-07 7.60807656e-07 7.60798628e-07 ... 7.60844159e-07
   7.60834965e-07 7.60825832e-07]
  ...
  [7.57482240e-07 7.57503419e-07 7.57524755e-07 ... 7.57423322e-07
   7.57441157e-07 7.57461227e-07]
  [7.57088410e-07 7.57111759e-07 7.57138485e-07 ... 7.57015632e-07
   7.57040577e-07 7.57065332e-07]
  [7.56607056e-07 7.56643601e-07 7.56680388e-07 ... 7.56530966e-07
   7.56556020e-07 7.56581372e-07]]
 [[5.84976528e-07 5.84976528e-07 5.84976528e-07 ... 5.84976528e-07
   5.84976528e-07 5.84976528e-07]
  [5.85080694e-07 5.85080831e-07 5.85080983e-07 ... 5.85080211e-07
   5.85080405e-07 5.85080533e-07]
  [5.85196801e-07 5.85197076e-07 5.85197379e-07 ... 5.85195844e-07
   5.85196124e-07 5.85196460e-07]
  ...
  [5.86448569e-07 5.86454287e-07 5.86460080e-07 ... 5.86431800e-07
   5.86437371e-07 5.86442908e-07]
  [5.86272733e-07 5.86279015e-07 5.86285173e-07 ... 5.86240342e-07
   5.86253994e-07 5.86266566e-07]
  [5.86108465e-07 5.86121208e-07 5.86134017e-07 ... 5.86047686e-07
   5.86067822e-07 5.86088111e-07]]]
```

### CAMS.getDataCAMS &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>

> [!NOTE]
> If a significant amount of data needs to be downloaded, we recommend to run CAMS.getDataCAMS [via Command Line Interface](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli).
> With this, you can download data sets in parallel - one set (_e.g._, one entire month) per Command Prompt.

```python
CAMS.getDataCAMS(dataType, years, months, days='All', hours=[0, 12], dataset=None, pressure_levels=[50, 100, 200, 400, 600, 800, 900, 1000], variables=None, bbox=[90, -180, -90, 180], mode=None, method='linear')
```
Download data from CAMS Global Greenhouse Gas Forecasts database.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ or _list of str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type(s) of data (e.g., 'mly' and/or 'dly').<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **years : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data (e.g., 2024).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **months : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month(s) of data (e.g., [4, 5]).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **days : _int_ or _list of int_ or _str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day(s) of month of data (e.g., [1, 5, 15]).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **hours : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Hour(s) of data (e.g., [0, 6, 12, 18]).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataset : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CAMS dataset name (e.g., "cams-global-greenhouse-gas-forecasts", "cams-global-ghg-reanalysis-egg4", "cams-global-atmospheric-composition-forecasts").<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pressure_levels : _int or List of int_, _optional, default: [50, 100, 200, 400, 600, 800, 900, 1000]_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A list of pressure levels to download (e.g., [100, 200, ..., 1000]).<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **variables : _str_ or _List of str_, _optional_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A list of variables to download (Allowed: "co", "co2", "ch4"). If None, uses the dataset defaults.<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _List_, _optional, default: [90, -180, -90, 180]_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earth's region of data, the bounding box.<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (upper_right_latitude, lower_left_longitude, lower_left_latitude, upper_right_longitude)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mode : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Mode for the download of data (Allowed: "add"). If "add", adds variable(s) to downloaded data. If None, downloads new data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **method : _str_, _optional_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method of interpolation (default: 'linear').<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **NPZ file in folder `data\CAMS\mly\` and/or `data\CAMS\dly\`**<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

<a name="CAMSMERRA2">**CAMS-MERRA2 atmospheric model**</a><br>
Combination of Copernicus Atmosphere Monitoring Service ([CAMS](#CAMS)) model and Modern-Era Retrospective analysis for Research and Applications Version 2 ([MERRA2](#MERRA2)).

### CAMSMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
instance_CAMSMERRA2 = CAMSMERRA2(dataType, y, m=None, d=None, bbox=(-180, -90, 180, 90), keys='All', phase='All', altArray=None, numAlt=50, surftrop=None, keysAsAttributes=False, showMessage=True)
```
Create an instance of `CAMSMERRA2` object.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'dly' - Daily data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly' - Monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'yly' - Yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cmly' - Combined monthly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'cyly' - Combined yearly data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ or _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _int_ for dataType 'dly', 'mly' and 'yly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y : _list of int_ [start_year, end_year] for dataType 'cmly' and 'cyly'.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _tuple_, _optional_, default: (-180, -90, 180, 90)_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Earth's region of data, the bounding box.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _str_ or _list of str_, _optional_, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Selection of phase of vertical profile composition.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'G' - Gas phase.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L' - Liquid phase.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-FW' - Liquid freshwater. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-SW' - Liquid seawater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'All' - All phases.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **altArray : _list of float_ or _np.ndarray of floats_, _optional_, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested altitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **numAlt : int_, _optional_, default: 50_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of altitude steps from 0.0 m to maximum tropopause altitude, if list of altitude is not given by the user with `altArray`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **surftrop : _str_ ('surface', 'tropopause'), _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Get concentration from 2-meters air following topography (`surftrop='surface'`) or tropopause height (`surftrop='tropopause'`).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keysAsAttributes : _bool_, _optional_, default: False_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set if keys from MERRA2 database are saved as object attributes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **showMessage : _bool_, _optional, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Boolean to set whether informative messages are displayed in Console.<p>
**Attributes:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.altitude : _list_ or _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric altitudes in meters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If surftrop = None, `.altitude` is a list of altitudes in meters from 0.0 to max tropopause altitude.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If surftrop = 'surfase', `.altitude` is the Earth's topography in meters. Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If surftrop = 'tropopause', `.altitude` is the tropopause altitude in meters. Shape: (latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.temperature : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric temperatures in Kelvin.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.pressure : _np.ndarray_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric pressures in Pascals.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Pi : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Partial pressure of compounds throughout the atmosphere in Pascals. `{'compound': [partial pressure]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_G : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gas concentration of compounds throughout the atmosphere in mol/L. `{'compound': [gas concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.MolPct_G : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gas molar percentage of compounds throughout the atmosphere in %mol. `{'compound': [molar percentage]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_LFW : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of compounds in freshwater aerosol throughout the atmosphere in mol/L. `{'compound': [liquid concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_LSW : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of compounds in seawater aerosol throughout the atmosphere in mol/L. `{'compound': [liquid concentration]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lat : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of latitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.lon : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of longitudes.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.compounds : _list_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atmospheric compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: (altitude, latitude, longitude).<p>

If _MERRA2_ data has already been downloaded, the user can create a new _CAMSMERRA2_ object with, at least, `dataType` and `y` arguments. The user can get the data calling the attributes defined above. Here is an example:
```python
from envdef import CAMSMERRA2

# Yearly data from 2020 was previously downloaded
newCAMSMERRA2 = CAMSMERRA2(dataType = 'yly', y = 2020, bbox = (-180, -90, -178.125, -88.5))
>>> print(newCAMSMERRA2.getAttributeNames())
['compounds', 'temperature', 'pressure', 'altitude', 'Pi', 'Ci_G', 'MolPct_G', 'Ci_LFW',
'Ci_LSW', 'lat', 'lon', 'environment', 'model']

newCAMSMERRA2 = CAMSMERRA2(dataType = 'yly', y = 2020, bbox = (-180, -90, -178.125, -88.5), keysAsAttributes = True)
>>> print(newCAMSMERRA2.getAttributeNames())
['compounds', 'temperature', 'pressure', 'altitude', 'Pi', 'Ci_G', 'MolPct_G', 'Ci_LFW',
'Ci_LSW', 'lat', 'lon', 'H', 'PS', 'PS_std', 'TROPPB', 'TROPPB_std', 'T2M', 'T2M_std',
'TROPT', 'TROPT_std', 'TROPH', 'TROPH_std', 'LR', 'LR_std', 'environment', 'model']

# Get data
>>> print(newCAMSMERRA2.lat)
[-90.  -89.5 -89.  -88.5]
>>> print(newCAMSMERRA2.lon)
[-180.    -179.375 -178.75  -178.125]
>>> print(newCAMSMERRA2.T2M)
[[224.31621 224.31621 224.31621 224.31621]
 [223.43001 223.42714 223.4242  223.42134]
 [223.60979 223.6157  223.62134 223.62642]
 [225.37038 225.3871  225.40935 225.43036]]
>>> print(newCAMSMERRA2.temperature)
# newCAMSMERRA2.temperature has NaN values because altitude < surface height ASL (topography).
[[[         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]]
 [[         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]
  [         nan          nan          nan          nan]]
...
[[208.89917863 208.89917863 208.89917863 208.89917863]
  [208.91566434 208.91529312 208.91491576 208.9145402 ]
  [209.0320008  209.03220517 209.0325211  209.03274575]
  [209.21578416 209.21842957 209.21895037 209.21945189]]
 [[208.4474782  208.4474782  208.4474782  208.4474782 ]
  [208.48575225 208.48551779 208.48527922 208.48503976]
  [208.59335294 208.59345043 208.59367159 208.59381527]
  [208.72229302 208.72450837 208.72443997 208.72438963]]]
>>> print(newCAMSMERRA2.Ci_LFW['O2'])
# newCAMSMERRA2.Ci_LFW['O2'] has NaN values because altitude < surface height ASL (topography).
[[[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
 [[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
...
 [[0.0005355  0.0005355  0.0005355  0.0005355 ]
  [0.00053454 0.00053452 0.00053449 0.00053447]
  [0.00053324 0.00053324 0.00053325 0.00053325]
  [0.00053258 0.00053257 0.00053256 0.00053254]]
 [[0.00052614 0.00052614 0.00052614 0.00052614]
  [0.00052481 0.00052478 0.00052475 0.00052472]
  [0.00052368 0.00052369 0.00052369 0.0005237 ]
  [0.00052402 0.00052402 0.00052402 0.00052401]]]
>>> print(newCAMSMERRA2.Ci_LFW['CH4'])
# newCAMSMERRA2.Ci_LFW['CH4'] has NaN values because altitude < surface height ASL (topography).
[[[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
 [[       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]
  [       nan        nan        nan        nan]]
...
 [[4.61909908e-09 4.61909908e-09 4.61909908e-09 4.61909908e-09]
  [4.61689825e-09 4.61696664e-09 4.61703609e-09 4.61710529e-09]
  [4.59911048e-09 4.59910210e-09 4.59907746e-09 4.59906631e-09]
  [4.57145217e-09 4.57108076e-09 4.57101823e-09 4.57095888e-09]]
 [[4.68624836e-09 4.68624836e-09 4.68624836e-09 4.68624836e-09]
  [4.68073763e-09 4.68078671e-09 4.68083658e-09 4.68088658e-09]
  [4.66392987e-09 4.66393717e-09 4.66392611e-09 4.66392672e-09]
  [4.64388328e-09 4.64356777e-09 4.64359123e-09 4.64361240e-09]]]
```

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

<a name="GWB">**General (or non-specific) Water Body (GWB)**</a><br>
Definition of a general (or non-specific) water body (abbreviated as GWB). A general (or non-specific) water body refers to any generic, unnamed, or unspecific natural or artificial collection of water, such as a lake, river, puddle, pool, pond, or stream.

### GWB &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
instance_GWB = GWB(Ct, T=[298.15], pH=[7.0], salinity=[0.0], fluidType='ideal', methods=None, showMessage=True)
```
Create an instance of `GWB` object.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Composition of water in mol/L. `{'compound': [concentration]}`<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _float_ or _list of floats_, _optional, default: 298.15_ **<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_ or _list of floats_, _optional, default: 7.0_ **<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **salinity : _float_ or _list of floats_, _optional, default: 0.0_ **<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Salinity of water body.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **fluidType : _str_ ('ideal' or 'non-ideal'), _optional, default: ideal_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of fluid (ideal or non-ideal).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **methods : _dict_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method for coefficient activity estimation `{'compounds': 'methods'}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'DH-ext'    - Debye-Hückel equation extended version.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'SS'        - Setschenow-Shumpe equation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **showMessage : _bool_, _optional, default: True_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Boolean to set whether informative messages are displayed in Console.<p>
**Attributes:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.temperature : _float_ or _list of floats_, _optional, default: 298.15_ **<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.pH : _float_ or _list of floats_, _optional, default: 7.0_ **<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.salinity : _float_ or _list of floats_, _optional, default: 0.0_ **<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Salinity of water body.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.Ci_L : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Composition of water in mol/L. `{'compound': [concentration]}`<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.fluidType : _str_ ('ideal' or 'non-ideal'), _optional, default: ideal_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of fluid (ideal or non-ideal).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *.*methods : _dict_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method for coefficient activity estimation `{'compounds': 'methods'}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'DH-ext'    - Debye-Hückel equation extended version.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'SS'        - Setschenow-Shumpe equation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **.compounds : _list of str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of considered compounds.<p>

The `GWB` does not need data to be download. The user can create a new _GWB_ object with, at least, `Ct` argument. The user can get the data calling the attributes defined above. Here is an example:
```python
from envdef import GWB
from thermodynamics import ThSA

new_GWB = GWB(Ct = {'CO2': [6.007e-4], 'Mg+2': [2.325e-2], 'Ca+2': [4.596e-3], 'Na+': [1.833e-1], 'K+': [4.032e-3],
                    'Cl-': [2.555e-1], 'Ba+2': [4.545e-9], 'Sr+2': [3.432e-6], 'SO4-2': [1.244e-2], 'NO3-': [5.432e-6],
                    'O2': [6.250e-6], 'H2': [7.213e-4], 'CH4': [3.472e-4], 'HCOOH': [4.000e-6], 'CH3COOH': [7.000e-5],
                    'NH2CH2COOH': [5.400e-9], 'H2S': [2.000e-12], 'NH3': [3.600e-7], 'NO2-': [1.000e-8], 'N2': [1.000e-3]},
              fluidType = 'ideal',
              methods = {'CO2': 'SS', 'Mg+2': 'DH-ext', 'Ca+2': 'DH-ext', 'Na+': 'DH-ext', 'K+': 'DH-ext',
                         'F-': 'DH-ext', 'Cl-': 'DH-ext', 'Ba+2': 'DH-ext', 'Sr+2': 'DH-ext', 'SO4-2': 'DH-ext',
                         'NO3-': 'DH-ext', 'O2': 'SS', 'H2': 'SS', 'CH4': 'SS', 'HCOOH': 'DH-ext', 'CH3COOH': 'DH-ext',
                         'NH2CH2COOH': 'DH-ext', 'H2S': 'SS', 'NH3': 'SS', 'NO2-': 'DH-ext', 'N2': 'SS'},
              T = np.arange(32.0 + 273.15, 46 + 273.15, 2.0),
              pH = [3.0, 4.0])

# Compute and show non-standard Gibbs free energy of formate generation (ForGen) and acetate Gen (AcGen)
new_GWB.getDGr('microprony', ['ForGen', 'AcGen'], ['H2', 'H2'])
>>> print(newGWB.DGr)
{'ForGen_pH:3.0': array([5.69677612, 5.74490896, 5.79296237, 5.84094139, 5.88885094, 5.93669582, 5.98448072]),
 'AcGen_pH:3.0': array([-20.38694508, -20.14563383, -19.90438742, -19.66320868, -19.42210047, -19.18106572, -18.94010739]),
 'ForGen_pH:4.0': array([0.00188635, 0.0106277 , 0.01950309, 0.02852364, 0.03770105, 0.04704761, 0.05657621]),
 'AcGen_pH:4.0': array([-20.49800783, -20.26430799, -20.03102153, -19.79815343, -19.56570758, -19.3336868 , -19.10209272])}
```

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)
#

<a name="create-new-environment">**How to create a new environment (class or subclass)**</a><br>
New environment classes (_Parent environment_) or subclasses (_child environment_) can be created in `envdef.py` module. Remember that a _child environment_ has all attributes and methods of _parent environments_ (_i.e._, `Environment01` and `Environment02` of example below). The _child environment_ can have its own attributes and methods in addition to those of the _parent environments_.

· For a new _parent environment_, follow the example below:
```python
class EnvironmentName:
    """
    Class information.
    """
    # Method to initialize a new instance (if initialization is necessary)
    def __init__(self, parameters):
        self.attribute = value
        self._privateAttribute = value

    # Instance method
    def method_01(self, paramters):
        """
        Method information.
        """
        # Code of method_01()

    # Another instance method
    def method_02(self, paramters):
        """
        Method information.
        """
        # Code of method_02()

    # Private instance method
    def _privateMethod(self, paramters):
        """
        Method information.
        """
        # Code of _privateMethod()
```

For more information about OOP in Python, click [here](https://realpython.com/python3-object-oriented-programming/). For more information of Constructors in Python (a special method that is called automatically when an object is created from a class), click [here](https://www.geeksforgeeks.org/constructors-in-python/).

· For a new _child environment_, follow the exemple below:
```python
class Environment01:
    """
    Class information.
    """
    # Method to initialize a new instance of Environment01 (if initialization is necessary)
    def __init__(self, parameters):
        self.attribute = value
        self._privateAttribute = value

    # Instance method of Environment01
    def methodName(self, paramters):
        """
        Method information.
        """
        # Code of methodName()

    # Private instance method
    def _privateMethod(self, paramters):
        """
        Method information.
        """
        # Code of _privateMethod()

class Environment02:
    """
    Class information.
    """
    # Method to initialize a new instance of Environment02 (if initialization is necessary)
    def __init__(self, parameters):
        self.attribute = value

    # Instance method of Environment02
    def methodName(self, paramters):
        """
        Method information.
        """
        # Code of methodName()

class ChildEnvironment(Environment01, Environment02):
    """
    Class information.
    """

    # Method to initialize a new instance of ChildEnvironments (required)
    def __init__(self, parameters01, parameters02):
        Environment01.__init__(self, parameters01)
        Environment02.__init__(self, parameters02)
        self.attribute = value

    # Instance method of ChildEnvironment
    def methodChild(self, paramters):
        """
        Method information.
        """
        # Code of methodChild()   
```

For more information of Multiple Inheritance in Python, click [here](https://www.geeksforgeeks.org/multiple-inheritance-in-python/).<br>
A common error in Python Inheritance is `TypeError: got multiple values for keyword argument`. You can fix this error following the example above and/or consulting [this website](https://www.geeksforgeeks.org/how-to-fix-python-multiple-inheritance-generates-typeerror-got-multiple-values-for-keyword-argument/).

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

#### <ins>Ecosystem Analysis (EcoA)</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

#### <ins>Thermodynamic State Analysis (ThSA)</ins>
With the <ins>Thermodynamic State Analysis module (ThSA)</ins>, the user can identify which rections are feasible (_i.e.,_ exergonic) considering the environmental conditions (nonstandard conditions). It is important to take into account that because a given reaction is exergonic under a particular environmental conditions (ΔG<sub>r</sub><0), it does not necessary mean that organisms will be able to catalyze it and, if they can, they have enough energy for growth and/or maintenance. Remember that Gibbs energy quantifies the tendency of a chemical reaction to proceed in a particular direction. To determine the actual energy that organisms can take from a particular transformation and if these can growth and/or satisfy maintenance requirements, the [Biological Thermodynamic State Analysis (BioThSA)](#bio-thermodynamic-state-analysis-biothsa) should be used. To calculate the nonstandard Gibbs free energy of reaction (ΔG<sub>r</sub>), the Gibbs-Helmholtz-Nernst relationship is used. The Gibbs-Helmholtz-Nernst relationship considers the influnce of temperature, pH and concentrations of substrates and products on ΔG<sub>r</sub>.

The main and auxiliary functions to perform the ThSA are located in `thermodynamics.py` module, and organized in three classes:
- `ThSA`. Class with the functions to perform the ThSA and export the results.
- `ThP`. Class with the functions to get thermodynamic parameters and calculate Gibbs energies and enthalpies of reactions.
- `ThEq`. Class with the functions to compute thermodynamic equilibriums (such as Henry solubility and pH speciation).

### ThSA.exportDeltaGr &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThSA.exportDeltaGr(modeExport, typeRxn, input_, phase, T, pH=7.0, S=None, Ct=1.0, specComp=False,
                   altitude=False, fluidType='ideal', molality=True, methods=None, solvent='H2O',
                   asm='stoich', warnings=False, printDG0r=False, printDH0r=False)
```
Compute the nonstandard Gibbs free energy of reaction (ΔG<sub>r</sub>) along the given conditions.<br> 
Resultant ΔG<sub>r</sub> is plotted in Spyder or written in an Excel file.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **modeExport : _str ('plot' or 'Excel')_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Way in which ΔG<sub>r</sub> is returned.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'plot' - ΔG<sub>r</sub> is plotted in Spyder.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'Excel' - ΔG<sub>r</sub> is written in an Excel file.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; What reaction database is used, matching with CSV name in `reactions\` folder.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **input_ : _str_ or _list of strs_ or _ndarray of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested compound(s) or reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str ('G' or 'L')_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Phase in which reaction(s) occurs.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _float_, _list of floats_, _ndarray of floats_, _optional, default: 298.15_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values [K].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _list of floats_, _ndarray of floats_, _optional, default: 7.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **S : _float_, _list of floats_, _ndarray of floats_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Salinity [ppt].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _dict_, _optional, default: 1.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total concentrations of compounds `{'compounds': [concentrations]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **specComp : _bool_, _str_, _list of strs_, _ndarray of strs_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Bool_ if `input_` are compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _str_, _list of strs_ or _ndarray of strs_ if `input_` are reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **altitude : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Altitude in _coordenate y_ of plots.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **fluidType : _str_ ('ideal' or 'non-ideal'), _optional, default: ideal_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of fluid (ideal or non-ideal).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **molality : _bool_, _optional, default: True_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Select if activity units are in molality (True) or molarity (False).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **methods : _dict_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method for coefficient activity estimation `{'compounds': 'methods'}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'DH-ext'    - Debye-Hückel equation extended version.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'SS'        - Setschenow-Shumpe equation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **solvent : _str_, _optional, default: H2O_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Solvent name.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **asm : _str ('asm')_, _optional, default: asm (stoichiometric concentrations)_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Assumption to calculate concentration of products not present in the environment.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **printDG0r : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print in console the values of standard Gibbs free energy of reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **printDH0r : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print in console the values of standard enthalpy of reactions.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Spyder plot** (if _modeExport='plot'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Excel file in `results\` folder** (if _modeExport='Excel'_)<br>

### ThSA.getDeltaGr &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThSA.getDeltaGr(typeRxn, input_, phase, specComp=False, T=298.15, pH=7.0, S=None, Ct=1.0,
                fluidType='ideal', molality=True, methods=None, solvent='H2O', asm='stoich',
                warnings=False, printDG0r=False, printDH0r=False)
```
Compute the nonstandard Gibbs free energy of reaction (ΔG<sub>r</sub>) along the given conditions.<br> 
Return a n-dimension array with ΔG<sub>r</sub> values.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; What reaction database is used, matching with CSV name in `reactions\` folder.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **input_ : _str_ or _list of strs_ or _ndarray of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested compound(s) or reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str ('G' or 'L')_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Phase in which reaction(s) occurs.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **specComp : _bool_, _str_, _list of strs_, _ndarray of strs_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Bool_ if `input_` are compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _str_, _list of strs_ or _ndarray of strs_ if `input_` are reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _float_, _list of floats_, _ndarray of floats_, _optional, default: 298.15_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values [K].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _optional, default: 7.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **S : _float_, _list of floats_, _ndarray of floats_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Salinity [ppt].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _dict_, _optional, default: 1.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total concentrations of compounds `{'compounds': [concentrations]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **fluidType : _str_ ('ideal' or 'non-ideal'), _optional, default: ideal_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of fluid (ideal or non-ideal).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **molality : _bool_, _optional, default: True_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Select if activity units are in molality (True) or molarity (False).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **methods : _dict_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method for coefficient activity estimation `{'compounds': 'methods'}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'DH-ext'    - Debye-Hückel equation extended version.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'SS'        - Setschenow-Shumpe equation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **solvent : _str_, _optional, default: H2O_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Solvent name.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **asm : _str ('asm')_, _optional, default: asm (stoichiometric concentrations)_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Assumption to calculate concentration of products not present in the environment.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **printDG0r : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print in console the values of standard Gibbs free energy of reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **printDH0r : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print in console the values of standard enthalpy of reactions.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **DGr : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Nonstandard Gibbs free energy of reaction.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(Z)x(Y)x(X)x(reactions)_.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For atmosphere: (Altitude)x(Latitude)x(Longitude)x(compounds).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user.<br>

### ThEq.plotpHSpeciation &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThEq.plotpHSpeciation(compounds, pH, temperature)
```
Plot pH (or ion) speciation of requested compounds.<p> 
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compounds : _str_, _list of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH :  _list of floats_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of requested pH values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature :  _float_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Temperature value for pH speciation.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Spyder plot** <br>

### ThEq.pHSpeciation &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThEq.pHSpeciation(iCompound, pH, t, Ct, rAllConc=False)
```
Compute pH (or ion) speciation of selected compounds.<br> 
Return a n-dimension array with concentrations of all chemical species.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **iCompound : _str_, _list of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested compound.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested pH value.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature : _float_, _list of floats_ or _ndarray of floats_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Temperature(s) for pH speciation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _list of float_, _ndarray of floats_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of total concentrations of compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rAllConc : _bool_, _optinal, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of total concentrations of compounds.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rSpec : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentrations of chemical species.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(Z)x(Y)x(X)_, (if _rAllConc = False_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(Z)x(Y)x(X)x(species)_, (if _rAllConc = True_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For atmosphere: (Altitude)x(Latitude)x(Longitude)x(species).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Where species: [B], [B<sup>-</sup>], [B<sup>-2</sup>], [B<sup>-3</sup>].<br>

### ThEq.solubilityHenry &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThEq.solubilityHenry(compounds, wType='FW', t=None)
```
Compute gas-liquid equilibrium based on Henry's law.<br> 
Return a n-dimension array with Henry's law solubility constant(s) and an array with the `compounds` indices of parameters that are available.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compounds : _str_, _list of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **wType : _str ('FW' or 'SW')_, _optional, default: FW_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Water type (phase).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'FW' - Freshwater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'SW' - Seawater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature : _float_ or _ndarray of floats_, _optional, default: None** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Hs : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Henry's law solubility constant(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(Z)x(Y)x(X)x(compounds)_, (if _temperature != empty_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For atmosphere: (Altitude)x(Latitude)x(Longitude)x(compounds).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(1)x(compounds)_, (if _temperature = empty_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **notNaN : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The `compounds` indices of parameters that are available.<br>

### ThP.getThP &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThP.getThP(typeParam, compounds, phase)
```
Load the thermodynamic parameters from local databases (`db\` folder).<br>
Return a n-dimension array with thermodynamic parameters and an array with the `compounds` indices of parameters that are available.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeParam : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of parameter database.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compounds : _str_, _list of strs_, or _ndarray of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Phase of parameter requested. Phase value depends on the database. See local databases in `db\Excels\` folder.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Param : _float_ or _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of requested parameters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **notNaN : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The `compounds` indices of parameters that are available.<br>

### ThP.getDeltaG0r &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThP.getDeltaG0r(deltaG0f, mRxn)
```
Compute standard Gibbs free energy of reaction from the standard Gibbs free energy of formation of substrate and products.<br>
Return a n-dimension array with the values of standard Gibbs free energy of reactions.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **deltaG0f : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of the standard Gibbs free energy of formation of compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **deltaG0r : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of the standard Gibbs free energy of reactions.<br>

### ThP.getDeltaH0r &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThP.getDeltaH0r(deltaH0f, mRxn)
```
Compute standard enthalpy of reaction from the standard enhalpy of formation of substrate and products.<br>
Return a n-dimension array with the values of standard enthalpy of reactions.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **deltaH0f : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of the standard enthalpy of formation of compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **deltaH0r : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of the standard enthalpy of reactions.<br>

### ThP.getKeq &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThP.getKeq(compounds, mRxn, t, phase)
```
Compute the equilibrium constant from the standard Gibbs free energy of reaction.<br>
Return a n-dimension array with the values of equilibrium constants.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compounds : _str_, _list of strs_, or _ndarray of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature : _float_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Temperature value.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Phase of parameter requested. Phase value depends on the database. See local databases in `db\Excels\` folder.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Keq : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(Z)x(Y)x(X)x(compounds)_, (if _temperature != empty_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For atmosphere: (Altitude)x(Latitude)x(Longitude)x(compounds).<br>

### ThP.ionicStrength &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThP.ionicStrengt(composition)
```
Compute ionic strength of solution.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **composition : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Composition of solution (`{'compounds': [concentrations]}`).<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **I : _float_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ionic strength of solution.<br>

### ThP.activity &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThP.activity(methods, composition, T=None, pH=None, salinity=None, molality=True, solvent='H2O', selComp=None)
```
Compute activities of compounds.<br>
Return a dictionary with activities of compound(s) in solution and their chemical species (if pH is given).<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **methods : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method for coefficient activity estimation `{'compounds': 'methods'}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'DH-ext'    - Debye-Hückel equation extended version.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'SS'        - Setschenow-Shumpe equation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **composition : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Composition of solution (`{'compounds': [concentrations]}`).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _float_, _list of floats_, _ndarray of floats_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values [K].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **salinity : _float_, _list of floats_, _ndarray of floats_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Salinity [ppt].<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **molality : _bool_, _optional, default: True_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Select if activity units are in molality (True) or molarity (False).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **solvent : _str_, _optional, default: H2O_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Solvent name.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **selComp : _str_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; To obtain the activity of a specific compound (must be in `methods` and `composition`).<p>

**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **activity : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Activity values (`{'compounds': [activity]}`)<br>

#

The functions to obtain the required information from reaction databases are found in `reaction.py` module. The main function is `Reactions.getRxn`.
> [!NOTE]
> If a same name is used for a compound and a reaction, use `getRxnByComp()` or `getRxnByName()` instead.

### Reactions.getRxn &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Reactions.getRxn(typeRxn, input_='All', warnings=False)
```
Load reaction information from local databases (`reactions\` folder).<br>
Return a list with involving compounds, a n-dimension array with stoichiometric matrix of reaction and information of reaction.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reaction database.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **input_ : _str_ or _list of strs_**, _optional, default: 'All'_ <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested compound(s) or reaction(s). <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If `input_='All'`, the user gets all reactions from `typeRxn.csv`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rComp : _list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of involving compounds of reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user. The function always returns the abbreviation.<br>

### Reactions.getRxnByComp &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Reactions.getRxnByComp(typeRxn, compounds, warnings = False)
```
Load reaction information from local databases (`reactions\` folder) based on given compound name(s).<br>
Return a list with involving compounds, a n-dimension array with stoichiometric matrix of reaction and information of reaction.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reaction database.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compounds : _str_ or _list of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested compound(s). <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rComp : _list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of involving compounds of reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user. The function always returns the abbreviation.<br>

### Reactions.getRxnByName &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Reactions.getRxnByName(typeRxn, nameRxn, warnings=False)
```
Load reaction information from local databases (`reactions\` folder) based on given reaction name(s).<br>
Return a list with involving compounds, a n-dimension array with stoichiometric matrix of reaction and information of reaction.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reaction database.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **nameRxn : _str_ or _list of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested reaction(s). <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rComp : _list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of involving compounds of reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **mRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user. The function always returns the abbreviation.<br>

#

**·** _<ins>Kinetics</ins>_

The he main and auxiliary functions to calculate kinetic rates are located in `reactions.py` module as well. These are organized in two classes:
- `KinP`. Class with the functions to get the kinetic parameters from local databases.
- `KinRates`. Class with the functions to calculate the kinetic rates (biotic and abiotic transformations).

### KinP.getKinP &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
KinP.getKinP(paramDB, params, reaction, sample='All', comp=None)
```
Load the kinetic parameters from local databases (`kinetics\` folder).<br>
Return a dictionary with the requested kinetic parameters and an array with the `sample` names (rows of `typeParam.csv` file).<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **paramDB : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of parameter database, matching with csv name.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **params : _str or list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of requested parameters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **reaction : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested reaction.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **sample : _str or list_, _optional, default: 'All'_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested samples (rows of `typeParam.csv`).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **comp : _str or list_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Compounds of parameters associated to compounds (e.g., Km).<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictR : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dictionary with requested parameters.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; dictR = {'Name parameter': [values]}.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **sampleNames : _str or list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of samples (_i.e._, rows of `typeParam.csv`).<br>

### KinRates.getRs &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
KinRates.getRs(typeKin, paramDB, reactions, Ct, sample = 'All', pH = None, T = None):
```
Compute reaction rates as a function of limiting substrate, inhibitors and temperatures based on the chosen kinetic equation (`typeKin`) and temperature correlation (`Tcorr`).<br> 
Return a n-dimension array with the calculated kinetic rates and an array with the `sample` names (rows of `typeParam.csv` file).<br>

> [!NOTE]
> If you want to calculate the influence of T for one rate, use single value arrays in `Ct` dictionary with same shape as temperature:
> ```python
> import numpy as np
> T = np.array([[273.15, 278.15, 283.15],
>              [288.15, 293.15, 298.15],
>              [303.15, 308.15, 313.15]])
> Ct = {'Compound A': 1.0 * np.ones(T.shape);
>      'Compound B': 2.0 * np.ones(T.shape);
>      'Compound C': 3.0 * np.ones(T.shape)}
> Rs, combNames, orderComb = typeKin = 'MM-Arrhenius', paramDB = ['MM_AtmMicr', 'ArrhCor_AtmMicr'],
>               reactions = 'Rxn1', Ct = Ct, sample = 'All', pH = 8.0, T = T)
> ```

**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeKin : _str_ ('MM' or 'MM-Arrhenihus')** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of kinetic equations.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'MM': Michaelis-Menten equation.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'MM-Arrhenihus': Michaelis-Menten-Arrhenius equation.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **paramDB : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of parameter database, matching with csv name in `kinetics\` folder (without '.csv').<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Rxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested reaction.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _dict of np.ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of substrates, products and/or inhibitors {'compound': np.array([concentrations])}.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **sample : _str or list_, _optional, default: 'All'_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested samples (rows of `paramDB.csv`).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float or list_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; pH value.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If pH value is given, ion speciation of compounds is computed.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _np.ndarray_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperatures. It must have the same size as concentrations of `Ct` dictionary.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Rs : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Resultant substrate uptake rates. {'reaction': {'comb_1': [rates], ..., 'comb_n': [rates]}} <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **combNames : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Combination of samples (rows of `typeParam.csv`). {'reaction: [combinations]'}<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **orderComb : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Order of samples (rows of `typeParam.csv`) in `combNames` (e.g., 'MM - Arrhenius').<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#
<!--
#### <ins>Bio-Thermodynamic State Analysis (BioThSA)</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#
-->
#### <ins>Ecosystem modelling</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

## :clipboard: Instructions to use EcoSysEM platform via Command Line Interface (CLI)
1. Download .zip code. Last version: `v0.1` **(Pre-release)**. [Download release](https://github.com/soundslikealloy/EcoSysEM/archive/refs/tags/v0.1.zip).
2. Extract files to a destination (Recommendation - Desktop).
3. Open **Anaconda Prompt or Terminal**.
4. Go to the **Code folder<sup>2</sup>** using `cd` command (more info about [Using Terminal](https://docs.anaconda.com/ae-notebooks/user-guide/basic-tasks/apps/use-terminal/?highlight=Using%20Terminal)).
    &#09;<br><sup><sup>2</sup>Code folder: folder with `ecosysem_cmd.py` file (Folder: `EcoSysEM\ecosysem`). </sup>
5. Execute one of the **EcoSysEM** blocks/functions using the following command lines:
```
python ecosysem_cmd.py _arg1 value1 _arg2 value2 _arg3 value3
```
Where `_arg#` are the arguments of the funtion and `value#` are the values of `_arg#`. Once executed the above command line, an `input()` line will request what function will be executed (see below). The user first gives all the arguments and corresponding values and then select the function.
```
> Available functions: getDataMERRA2, getDataCAMS
>> Enter the function:
```
### Arguments list:
#### <ins>getDataMERRA2</ins>
<table border="0">
   <tr><td> _h<br>__help </b></td><td> Show help message and optional arguments.</b></td></tr>
   <tr><td> _dataType </td><td> [str or list] Type of data (dly: daily; mly: monthly; cmly: combined monthly; All: dly mly cmly).</td></tr>
   <tr><td> _y </td><td> [int] Year of requested data.</td></tr>
   <tr><td> _m </td><td> [int or list] Month(s) of requested data.</td></tr>
   <tr><td> _d </td><td> [str or int] (Default: 'All') Last day of month of requested data. With 'All' get the whole month.</td></tr>
   <tr><td> _product </td><td> [str] (Default: 'M2I1NXASM') Product of data (section of MERRA2 database).</td></tr>
   <tr><td> _version </td><td> [str] (Default: '5.12.4') Version of data.</td></tr>
   <tr><td> _bbox </td><td> [tuple] (Default: '-180 -90 180 90) Earths region of data, the bounding box `-bbox lower_left_lon lower_left_lat upper_right_lon upper_right_lat`.</td></tr>
   <tr><td> _var </td><td> [list of str] (Default: PS T2M TROPT TROPPB) List of requested variables.</td></tr>
</table>

List and tuples are given without `[]` or `()`, and elements are separated by space. Strings are given without `' '` or `" "`. 
For example: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_dataType cmly mly` => `dataType = ['cmly' 'mly']` <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_y 2024 2025` => `year = [2024 2025]` <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_var PS T2M TROPT TROPPB` => `var = ['PS', 'T2M', 'TROPT', 'TROPPB']` <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_bbox -180 -90 -178.125 -86.5` => `bbox = (-180, -90, -178.125, -86.5)` <br>

More examples below:
```
python ecosysem_cmd.py _dataType mly _y 2021 _m 4
python ecosysem_cmd.py _dataType mly _y 2021 _m 4 _var PS TROPPB T2M TROPT H TROPH LR
python ecosysem_cmd.py _dataType dly mly _y 2021 _m 4 _var PS TROPPB T2M TROPT H TROPH LR
python ecosysem_cmd.py _dataType All _y 2021 2022 2023 _m 1 4 7 10 _bbox -180 -90 -178.125 -86.5
```

#### <ins>getDataCAMS</ins>
<table border="0">
   <tr><td> _h<br>__help </b></td><td> Show help message and optional arguments.</b></td></tr>
   <tr><td> _dataType </td><td> [str] Type(s) of data ('mly', 'dly').</td></tr>
   <tr><td> _y </td><td> [int or list] Year(s) of requested data.</td></tr>
   <tr><td> _m </td><td> [int or list] Month(s) of requested data.</td></tr>
   <tr><td> _d </td><td> [int or list or str 'All'] (Default: 'All') Day(s) of month of requested data. With 'All' get the whole month.</td></tr>
   <tr><td> _hr </td><td> [int or list of int] (Default: '0 12') Hour(s) of requested data.</td></tr>
   <tr><td> _dset </td><td> [str] Name of dataset ('cams-global-greenhouse-gas-forecasts', 'cams-global-ghg-reanalysis-egg4', 'cams-global-atmospheric-composition-forecasts').</td></tr>
   <tr><td> _pressure </td><td> [int or list] (Default: '[50, 100, 200, 400, 600, 800, 900, 1000]') Pressure levels to download.</td></tr>
   <tr><td> _bbox </td><td> [list] (Default: '90 -180 -90 180') Earth's region of data, the bounding box `-bbox upper_right_lat lower_left_lon lower_left_lat upper_right_lon`.</td></tr>
   <tr><td> _mode </td><td> [str] Mode of download ('add').</td></tr> 
   <tr><td> _method </td><td> [str] (Default: 'linear') Method of interpolation in 'add' mode.</td></tr> 
</table>

List and tuples are given without `[]` or `()`, and elements are separated by space. Strings are given without `' '` or `" "`. 
For example: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_y 2024 2025` => `year = [2024 2025]` <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_dset cams-global-greenhouse-gas-forecasts` => `dataset = 'cams-global-greenhouse-gas-forecasts'` <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_pressure 200 400` => `pressure_levels = [200, 400]` <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`_bbox 90 -180 -90 180` => `bbox = [90, -180, -90, 180]` <br>

More examples below:
```
python ecosysem_cmd.py _dataType mly _y 2024 _m 4
python ecosysem_cmd.py _dataType mly _y 2024 _m 4 _d 1 15
python ecosysem_cmd.py _dataType mly dly _y 2024 _m 4 5 _d 1 15 _pressure 200 300
python ecosysem_cmd.py _dataType mly _y 2024 _m 4 5 6 7 8 _bbox 90 -180 -90 180
```

[🔼 Back to **Contents**](#readme-contents)

## Function Navigation
#### · <ins>General functions for all environmental models (Environment object)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.loadData](#environmentloaddata---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.combData](#environmentcombdata---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.getAttributeNames](#environmentgetattributenames---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.getDGr](#getdgr---back-to-function-navigation)<br>

#### · <ins>ISA (Atmosphere)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA](#isa---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.setComposition](#isasetcomposition---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.plotTandP](#isaplottandp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.plotCompsProfilesISA](#isaplotcompsprofilesisa---back-to-function-navigation)<br>

#### · <ins>MERRA2 (Atmosphere)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2](#merra2---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.getDataMERRA2](#merra2getdatamerra2---back-to-function-navigation)<br>

#### · <ins>ISAMERRA2 (Atmosphere)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISAMERRA2](#isamerra2---back-to-function-navigation)<br>

#### · <ins>CAMS (Atmosphere)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [CAMS](#cams---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [CAMS.getDataCAMS](#camsgetdatacams---back-to-function-navigation)<br>

#### · <ins>CAMSMERRA2 (Atmosphere)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [CAMSMERRA2](#camsmerra2---back-to-function-navigation)<br>

#### · <ins>GWB (Hydrosphere)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [GWB](#GWB---back-to-function-navigation)<br>

#### · <ins>Thermodynamic equilibrium (ThEq)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThEq.plotpHSpeciation](#theqplotphspeciation---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThEq.pHSpeciation](#theqphspeciation---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThEq.solubilityHenry](#theqsolubilityhenry---back-to-function-navigation)<br>

#### · <ins>Thermodynamic parameters (ThP)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getThP](#thpgetthp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getDeltaG0r](#thpgetdeltag0r---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getDeltaH0r](#thpgetdeltah0r---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getKeq](#thpgetkeq---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.ionicStrength](#thpionicstrength---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.activity](#thpactivity---back-to-function-navigation)<br>

#### · <ins>Thermodynamic State Analysis (ThSA)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThSA.exportDeltaGr](#thsaexportdeltagr---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThSA.getDeltaGr](#thsagetdeltagr---back-to-function-navigation)<br>

#### · <ins>Reactions (Reactions)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Reactions.getRxn](#reactionsgetrxn---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Reactions.getRxnByComp](#reactionsgetrxnbycomp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Reactions.getRxnByName](#reactionsgetrxnbyname---back-to-function-navigation)<br>

#### · <ins>Kinetic parameters (KinP)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [KinP.getKinP](#kinpgetkinp---back-to-function-navigation)<br>

#### · <ins>Kinetic rates (KinRates)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [KinRates.getRs](#kinratesgetrs---back-to-function-navigation)<br>

[🔼 Back to **Contents**](#readme-contents)

## Error List
#### <ins>MERRA2.getDataMERRA2</ins>
**> _Invalid credentials (Earthacces username and/or password)_**<br>
```
[...]
earthaccess.exceptions.LoginAttemptFailure: Authentication with Earthdata Login failed with:
{"error":"invalid_credentials","error_description":"Invalid user credentials"}
```
**· Solution**: introduce valid credentials.<pr>

**> _Server Error (MERRA-2 server)_**<br>
```
[...]
RuntimeError: {"errors":["An Internal Error has occurred."]}
```
**· Solution 1**: accept all _end-user license aggrements_ (EULAs), if not yet done. You can find them after [login](https://urs.earthdata.nasa.gov/home) with your earthaccess account in <ins>EULAs</ins> -> <ins>Accept New EULAs</ins>.<br>
**· Solution 2**: restart the Anaconda prompt (open a new Anaconda prompt) or Spyder console (Ctrl + D) and run the code again.

#### <ins>CAMS.getDataCAMS</ins>
**> _Invalid temporary files_**<br>
```
[...]
ValueError: Failed for yyyy-mm-dd: Error occurred: CAMS_yyyy_m_d.nc, Failed to decode variable 'forecast_period': failed to prevent overwriting existing key dtype in attrs on variable 'forecast_period'. This is probably an encoding field used by xarray to describe how a variable is serialized. To proceed, remove this key from the variable's attributes manually.
```
**· Solution**: delete `/temporary` folder in `data/CAMS/` folder.<pr>

[🔼 Back to **Contents**](#readme-contents)

__________________________________________________

## Contact

**Eloi Martinez-Rabert**. :envelope: eloi.mrp@gmail.com

[🔼 Back to **Top**](#ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### References

[^1]: International Organization for Standardization, Standard Atmosphere, ISO 2533:1975, 1975. [Link to ISO](https://www.iso.org/standard/7472.html).
[^2]: The Atmosphere. Introduction to the Atmosphere. From National Oceanic and Atmospheric Administration (NOAA). Last updated: 2 July, 2024 [Link to NOAA](https://www.noaa.gov/jetstream/atmosphere).
[^3]: Goody, R. M., & Yung, Y. L. (1989). Atmospheric Radiation: Theoretical Basis. Oxford University Press. [https://doi.org/10.1093/oso/9780195051346.001.0001](https://doi.org/10.1093/oso/9780195051346.001.0001).
[^4]: [Acker et al. (2005)](https://www.sciencedirect.com/science/article/pii/S0169809504001462); [Hanke et al. (2003)](https://acp.copernicus.org/articles/3/417/2003/); [Lu et al. (2019)](https://acp.copernicus.org/articles/19/1971/2019/); [Zhong et al. (2023)](https://acp.copernicus.org/articles/23/14761/2023/)
