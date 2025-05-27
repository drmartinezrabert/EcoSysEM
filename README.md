# EcoSysEM platform
Eco-System Evaluation &amp; Modelling

*· Contributors: Eloi Martinez-Rabert* 

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
        - Ecosystem Analysis (EcoA) 🚧 | [GO](#ecosystem-analysis-ecoa)
        - Thermodynamic State Analysis (ThSA) | [GO](#thermodynamic-state-analysis-thsa)
        - Bio-Thermodynamic State Analysis (BioThSA) 🚧 | [GO](#bio-thermodynamic-state-analysis-biothsa)
        - Ecosystem modelling 🚧 | [GO](#ecosystem-modelling)
-  Instructions to use EcoSysEM platform via Command Line Interface (CLI) 🚧 | [GO](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli)
-  Function Navigation | [GO](#function-navigation)
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
[🔼 Back to **Contents**](#readme-contents)

____________________________

## :clipboard: Instructions for downloading and setting up EcoSysEM platform
1. Download .zip code. Last version: `v0.1` **$${\color{orange}\textbf{(Pre-release)}}$$**. [Download release](https://github.com/soundslikealloy/EcoSysEM/archive/refs/tags/v0.1.zip).
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
  ├── envdef.py 
  │      ├── ISA
  │      │    ├── .ISAaltitude
  │      │    ├── .ISAtemperature
  │      │    ├── .ISApressure
  │      │    ├── .compounds
  │      │    ├── .compositions
  │      │    ├── setComposition
  │      │    ├── selectAltitude
  │      │    ├── dictConcISA
  │      │    ├── plotTandPISA
  │      │    └── plotCompsProfilesISA
  │      ├── MERRA2
  │      │    ├── getDataMERRA2
  │      │    ├── combDataMERRA2
  │      │    ├── dictMERRA2
  │      │    ├── keysMERRA2
  │      │    └── deleteKeyMERRA2
  │      └── ISAMERRA2 {subclass of ISA & MERRA2}
  ├── reactions.py
  │      ├── KinP
  │      │    └── getKinP
  │      ├── KinRates
  │      │    ├── getRs
  │      │    ├── rMM
  │      │    └── arrhCorr
  │      └── Reactions
  │           ├── getRxn
  │           ├── getRxnByComp
  │           └── getRxnByName
  └── thermodynamics.py 
         ├── ThP
         │    ├── getThP
         │    ├── getDeltaG0r
         │    ├── getDeltaH0r
         │    └── getKeq
         ├── ThEq
         │    ├── solubilityHenry
         │    ├── pHSpeciation
         │    └── plotpHSpeciation
         └── ThSA
              ├── getDeltaGr
              └── exportDeltaGr
```

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### Fundamentals and usage
This section clarifies concepts, design decisions and technical details of this package. **EcoSystem platform** is constituted by four main units:
- Environment definition and instance calling | [GO](#environment-definition-and-instance-calling)
  - Ideal Earth's atmosphere (International Standard Atmosphere, ISA) | [GO](#ISA)
  - Modern-Era Retrospective analysis for Research and Applications, Version 2 (MERRA-2) | [GO](#MERRA2)
  - ISA-MERRA2 atmospheric model | [GO](#ISAMERRA2)
  - How to create a new environment (class or subclass) | [GO](#create-new-environment)
- Ecosystem Analysis (EcoA) 🚧 | [GO](#ecosystem-analysis-ecoa)
- Thermodynamic State Analysis (ThSA) | [GO](#thermodynamic-state-analysis-thsa)
- Bio-Thermodynamic State Analysis (BioThSA) 🚧 | [GO](#bio-thermodynamic-state-analysis-biothsa)
- Ecosystem modelling 🚧 | [GO](#ecosystem-modelling)

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Environment definition and instance calling</ins>
One of the advantages of Python is that supports both **Object-Oriented Programming (OOP)** and functional programming paradigms. The definition of environments are based on <ins>OOP paradigm</ins>. OOP is based on the following four principles: _Encapsulation_, _Inheritance_, _Abstraction_ and _Polymorphism_. For more info about OOP principles, click [here](https://fluxtech.me/blog/object-oriented-programming-vs-functional-programming/). 
- _Encapsulation_ principle allows to hide he internal state and behaviour of an object, and the object can only be acessed through a well-defined interface. With this, the user can change the properties of the environment without affecting the code hat uses the object.
- _Inheritance_ principle allows a new class to be defined based on an existing class, inheriting its attributes and methods.
- _Abstraction_ principle makes possible to work with objecs of a class without knowing the details of their implementation, which can make the code more robust and less error-prone.
- _Polymorphism_ principle enables the use of a common interface for different classes, making it possible to write code that can work with object of different types without knowing their specific class.

The benefits of OOP are _i_) organization, _ii_) state definition and tracking, _iii_) encapsulation of proceudre and data (_i.e.,_ specific functions and data can be stored together in a single class), _iv_) inheritance (making development more efficient and easier to maintain). For more information about OOP in Python, click [here](https://realpython.com/python3-object-oriented-programming/).

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

To create a new _ISA_ object (_i.e.,_ instantiate the class `ISA`), the instance attributes `layers`, `H2O`, `pH` and `resolution` are necessary:
- `layers`. Selection of atmosphere layers defined by ISA model[^1]. This attribute can be 'All' (_string_), an _integer_ from 0 to 7 or a _list_ of integers.
- `H2O`. Water content of atmosphere. This attribute must be a _float_ from 0.0 to 0.04.
- `pH`. pH of atmosphere. This attribute must be a _float_.
- `resolution`. Resolution of altitude array, that is, the size of altitude nodes per layer (in m). This attribute must be an _integer_.

Once a new _ISA_ object is created, a specific region of ISA can be selected using `ISA.selectRegion`. Compound concentrations (`ISA.getDictConc` or `ISA.getVerticalProfiles`) will be calculated using the new region defined. Here is an example:
```python
from envdef import ISA
import numpy as np

newISA = ISA(0, 0.00, None, 500)

>>> print(newISA.ISAaltitude)
[    0.   500.  1000.  1500.  2000.  2500.  3000.  3500.  4000.  4500.
  5000.  5500.  6000.  6500.  7000.  7500.  8000.  8500.  9000.  9500.
 10000. 10500. 11000.]
>>> print(newISA.ISAtemperature)
[ 15.    11.75   8.5    5.25   2.    -1.25  -4.5   -7.75 -11.   -14.25
 -17.5  -20.75 -24.   -27.25 -30.5  -33.75 -37.   -40.25 -43.5  -46.75
 -50.   -53.25 -56.5 ]

newISA.selectRegion([1000, 3500])

>>> print(newISA.ISAaltitude)
[1000. 1500. 2000. 2500. 3000. 3500.]
>>> print(newISA.ISAtemperature)
[ 8.5   5.25  2.   -1.25 -4.5  -7.75]

C = newISA.getDictConc('L-SW')

# np.round(a, decimals) -> NumPy function: Evenly round to the given number of decimals.
>>> print(np.round(C['O2'], 7)) 
[0.0003134, 0.0003138, 0.0003144, 0.0003153, 0.0003164, 0.0003178]
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

### ISA.selectAltitude &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.selectAltitude(selAlt)
```
Modify `.altitude`, `.temperature` and `.pressure` of `ISA` subclass based on the minimum and maximum altitude given.<p>
**Parameters:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **selAlt : _int, float or list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Minimum and maximum altitude of the new atmosphere region.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If _selAlt_ is a _int_ or _float_: minAlt = 0; maxAlt = selAlt.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If _selAlt_ is a _list_: [minAlt, maxAlt].<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **None**<br>

### ISA.dictConcISA &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.dictConcISA(phase, compound=None)
```
Return vertical profiles in _format=dict_ of selected compounds or all compounds of `ISA` subclass.<p>
**Parameters:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str ('G', 'L-FW', 'L-SW', 'L' or 'All')_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The desired phase of compound concentration.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'G' - Gas.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-FW' - Liquid freshwater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L-SW' - Liquid seawater.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L' - Both liquid phases (L-FW, L-SW).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'All' - All phases (G, L-FW, L-SW).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compound : _str_ or _list of strs_, _optional, default: None_ (all compounds)**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The desired compound(s) of `ISA` subclass.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictPi, dictCi_G : _dict_** (if _phase='G'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictCi_LFW : _dict_** (if _phase='L-FW'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictCi_LSW : _dict_** (if _phase='L-SW'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictCi_LFW, dictCi_LSW : _dict_** (if _phase='L'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictPi, dictCi_G, dictCi_LFW, dictCi_LSW : _dict_** (if _phase='All'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; dictPi : Parcial pressure of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; dictCi_G : Concentration in gas of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; dictCi_LFW : Concentration in liquid (freshwater) of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; dictCi_LSW : Concentration in liquid (seawater) of desired compounds.<br>

### ISA.plotTandP_ISA &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.plotTandP_ISA()
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

To create a new _MERRA2_ object (_i.e.,_ instantiate the class `MERRA2`), no instance attributes are necessary. Once a new _MERRA_ object is created, the available data can be downloaded and combined using `MERRA2.getDataMERR2` and `MERRA2.combDataMERRA2`, respectively. The data wil be saved in the folder `data\MERRA2\`. The data is saved in .npz file format (more info [here](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html). Once the data is downladed, the user can obtain the data with `MERRA2.dictMERRA2` function, see the parameters of data with `MERRA2.keysMERRA2`, or delete existing keys with `MERRA2.deleteKeyMERRA2`. Here is an example:
```python
from envdef import MERRA2

newMERRA2 = MERRA2()

# Get monthly data from online databases
## (Default arguments: product = 'M2I1NXASM', version = '5.12.4', var = ['PS', 'T2M', 'TROPT', 'TROPPB'])
newMERRA2.getDataMERRA2(years = 1995, months = 1, days = 'All', bbox = (-180, -90, -178.125, -88.5))

# See keys (_e.i._, variables names) of downloaded data
keys = newMERRA2.keysMERRA2(dataType = 'mly', y = 1995, m = 1)

>>> print(key)
['lat', 'lon', 'PS', 'PS_std', 'T2M', 'T2M_std', 'TROPT', 'TROPT_std',
'TROPPB', 'TROPPB_std', 'H', 'H_std', 'TROPH', 'TROPH_std', 'LR', 'LR_std']

# See data
data = newMERRA2.dictMERRA2(dataType = 'mly', y = 1995, m = 1, keys = ['lat', 'lon', 'T2M'])

>>> print(data)
{'lat': array([-90. , -89.5, -89. , -88.5]),
'lon': array([-180.   , -179.375, -178.75 , -178.125]),
'T2M': array([[244.29813, 244.29813, 244.29813, 244.29813],
              [243.8813 , 243.88293, 243.88455, 243.88733],
              [244.05316, 244.0517 , 244.04942, 244.04617],
              [245.6959 , 245.69054, 245.68338, 245.67622]], dtype=float32)}
>>> print(data['T2M'])
[[244.29813 244.29813 244.29813 244.29813]
 [243.8813  243.88293 243.88455 243.88733]
 [244.05316 244.0517  244.04942 244.04617]
 [245.6959  245.69054 245.68338 245.67622]]
```

### MERRA2.getDataMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>

> [!NOTE]
> If a significant amount of data needs to be downloaded, we recommend to run MERRA2.getDataMERRA2 [via Command Line Interface](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli).
> With this, you can download data sets in parallel - one set (_e.g._, one entire month) per Command Prompt.

```python
MERRA2.getDataMERRA2(years, months, days='All', product='M2I1NXASM', version='5.12.4', bbox=(-180, -90, 180, 90), var=['PS', 'T2M', 'TROPT', 'TROPPB'], daily=False)
```
Download data from MERRA2 database.<p>
**Parameters:**<br>
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **daily : _bool_, _optional, default: False_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Save daily data or not.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **NPZ file in folder `data\MERRA2\mly\` and `data\MERRA\dly\` if _daily = True_**<br>

### MERRA2.combDataMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
MERRA2.combDataMERRA2(years, month, dataType, keys = 'All')
```
Get average and standard deviation from a group of data.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **years : _list of int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Years of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **month : _int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ ('mly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **NPZ file in folder `data\MERRA2\cmly\`**<br>

### MERRA2.selectRegion &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
MERRA2.selectRegion(data, bbox)
```
Get average and standard deviation from a group of data.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **data : _dict_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data in dictionary form.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **bbox : _tuple_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested coordinates, the bounding box.<br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataSel : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data from requested region.<br>

### MERRA2.dictMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
MERRA2.dictMERRA2(dataType, y, m, d=None, keys='All')
```
Get data in dictionary form.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ ('mly', 'cmly', 'dly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly': monthly; 'cmly': combined monthly.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ ('mly' and 'dly') or _list of int_ ('cmly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_, _optional, default: 'All'_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of requested variables.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dictVar : _dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dictionary with requested variables.<br>

### MERRA2.keysMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
MERRA2.keysMERRA2(dataType, y, m, d=None)
```
Get variable list of data.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ ('mly', 'cmly', 'dly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly': monthly; 'cmly': combined monthly.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ ('mly' and 'dly') or _list of int_ ('cmly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Variable list of data.<br>

### MERRA2.deleteKeyMERRA2 &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
MERRA2.deleteKeyMERRA2(keys, dataType, y, m, d=None)
```
Delete variable(s) from data.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **keys : _list of str_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of variables to be deleted.<p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **dataType : _str_ ('mly', 'cmly', 'dly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'mly': monthly; 'cmly': combined monthly.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **y : _int_ ('mly' and 'dly') or _list of int_ ('cmly')**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Year(s) of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **m : _int_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Month of data.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **d : _int_, _optional, default: None_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Day of data.<p>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

New inhereted classes (also known as _subclasses_) can be created inhereting the attributes from existing classes (_e.g.,_ ISA or MERRA2). New attributes and function can be defined in these sublcasses, which will belong only to the subclass in question.

<a name="ISAMERRA2">**ISA-MERRA2 atmospheric model**</a><br>
Combination of International Standard Atmosphere ([ISA](#ISA)) model and Modern-Era Retrospective analysis for Research and Applications Version 2 ([MERRA2](#MERRA2)).

Because `ISAMERRA2` sublcass is a multiple-inhereted class of `ISA` and `MERRA2` classes, this has all attributes and methods of parent classes (_i.e._, `ISA` and `MERRA2`). All functions of `ISA` and `MERRA2` classes are summarized in [EcoSysEM package layout](#ecosysem-package-layout). `ISAMERRA2` has no new attributes or functions. 

To create a new _ISAMERRA2_ object (_i.e.,_ instantiate the class `ISAMERRA2`), no instance attributes are necessary. The attributes of `ISA` class are by default: layers = 0, H2O = 0.0, pH = 8.0, resolution = 1000. These attributes can be modified when the instance of `ISAMERRA2` is created. Here is an example:
```python
from envdef import ISAMERRA2

newISAMERRA2 = ISAMERRA2(layers=[0], resolution = 200)

# Get temperature profile from ISA model
>>> print(newISAMERRA2.ISAtemperature)
[ 15.   13.7  12.4  11.1   9.8   8.5   7.2   5.9   4.6   3.3   2.    0.7
  -0.6  -1.9  -3.2  -4.5  -5.8  -7.1  -8.4  -9.7 -11.  -12.3 -13.6 -14.9
 -16.2 -17.5 -18.8 -20.1 -21.4 -22.7 -24.  -25.3 -26.6 -27.9 -29.2 -30.5
 -31.8 -33.1 -34.4 -35.7 -37.  -38.3 -39.6 -40.9 -42.2 -43.5 -44.8 -46.1
 -47.4 -48.7 -50.  -51.3 -52.6 -53.9 -55.2 -56.5]

# Get atmospheric composition from ISA model
>>> print(newISAMERRA2.compositions)
{'N2': 0.78084, 'O2': 0.20946, 'Ar': 0.00934, 'CO2': 0.000426, 'Ne': 1.8182e-05,
'He': 5.24e-06, 'CH4': 1.92e-06, 'Kr': 1.14e-06, 'H2': 5.5e-07, 'N2O': 3.3e-07,
'CO': 1e-07, 'Xe': 9e-08, 'O3': 7e-08, 'NO2': 2e-08, 'SO2': 1.5e-08, 'I2': 1e-08,
'NH3': 6e-09, 'HNO2': 1e-09, 'HNO3': 1e-09, 'H2S': 3.3e-10}

# See keys (_e.i._, variables names) of downloaded data
keys = newISAMERRA2.keysMERRA2(dataType = 'mly', y = 1995, m = 1)

>>> print(key)
['lat', 'lon', 'PS', 'PS_std', 'T2M', 'T2M_std', 'TROPT', 'TROPT_std',
'TROPPB', 'TROPPB_std', 'H', 'H_std', 'TROPH', 'TROPH_std', 'LR', 'LR_std']

# See data
data = newISAMERRA2.dictMERRA2(dataType = 'mly', y = 1995, m = 1, keys = ['lat', 'lon', 'T2M'])

>>> print(data)
{'lat': array([-90. , -89.5, -89. , -88.5]),
'lon': array([-180.   , -179.375, -178.75 , -178.125]),
'T2M': array([[244.29813, 244.29813, 244.29813, 244.29813],
              [243.8813 , 243.88293, 243.88455, 243.88733],
              [244.05316, 244.0517 , 244.04942, 244.04617],
              [245.6959 , 245.69054, 245.68338, 245.67622]], dtype=float32)}
>>> print(data['T2M'])
[[244.29813 244.29813 244.29813 244.29813]
 [243.8813  243.88293 243.88455 243.88733]
 [244.05316 244.0517  244.04942 244.04617]
 [245.6959  245.69054 245.68338 245.67622]]
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
ThSA.exportDeltaGr(modeExport, T, pH, phase, typeRxn, input_, specComp=False, Ct=1.0, Ct_associated=None, asm='stoich', warnings=False)
```
Compute the nonstandard Gibbs free energy of reaction (ΔG<sub>r</sub>) along the given conditions.<br> 
Resultant ΔG<sub>r</sub> is plotted in Spyder or written in an Excel file.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **modeExport : _str ('plot' or 'Excel')_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Way in which ΔG<sub>r</sub> is returned.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'plot' - ΔG<sub>r</sub> is plotted in Spyder.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'Excel' - ΔG<sub>r</sub> is written in an Excel file.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _float_, _list of floats_, _ndarray of floats_ or _None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values. _None_ if standard temperature.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _list of floats_, _ndarray of floats_ or _None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH values. _None_ if standard pH.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str ('G' or 'L')_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Phase in which reaction(s) occurs.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'G' - Gas.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 'L' - Liquid.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeRxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; What reaction database is used, matching with CSV name in `reactions\` folder.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **input_ : _str_ or _list of strs_ or _ndarray of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of requested compound(s) or reaction(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **specComp : _bool_, _str_, _list of strs_, _ndarray of strs_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _Bool_ if `input_` are compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _str_, _list of strs_ or _ndarray of strs_ if `input_` are reactions.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _dict_, _optional, default: 1.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total concentrations of compounds `{'compounds': [concentrations]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct_associated : _str ('x', 'y' or 'xy')_ or _None_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set if concentrations are associated with temperature ('y'), pH ('x') or both.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **asm : _str ('asm')_, _optional, default: asm (stoichiometric concentrations)_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Assumption to calculate concentration of products not present in the environment.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Spyder plot** (if _modeExport='plot'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Excel file in `results\` folder** (if _modeExport='Excel'_)<br>

### ThSA.getDeltaGr &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThSA.getDeltaGr(typeRxn, input_, phase, specComp=False, T=298.15, Ct=1.0, pH=7.0, asm='stoich', warnings=False)
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _dict_, _optional, default: 1.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total concentrations of compounds `{'compounds': [concentrations]}`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _list of floats_, _ndarray of floats_, _optional, default: 7.0_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of pH values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **asm : _str ('asm')_, _optional, default: asm (stoichiometric concentrations)_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Assumption to calculate concentration of products not present in the environment.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **warnings : _bool_, _optional, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Display function warnings.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **DGr : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Nonstandard Gibbs free energy of reaction.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(reactions)x(temperature)x(pH)x(total concentration)_.<br>
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
ThEq.pHSpeciation(compounds, pH, temperature, Ct, rAllConc=False)
```
Compute pH (or ion) speciation of selected compounds.<br> 
Return a n-dimension array with concentrations of all chemical species.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compounds : _str_, _list of strs_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pH : _float_, _list of floats_ or _ndarray of floats_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of requested pH values.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature : _float_, _list of floats_ or _ndarray of floats_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Temperature(s) for pH speciation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _list of float_, _ndarray of floats_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of total concentrations of compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rAllConc : _bool_, _optinal, default: False_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of total concentrations of compounds.<p>
**Returns:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rSpec : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentrations of chemical species.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(pH)x(total concentration)x(temperature)x(compounds)_, (if _rAllConc = False_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(species)x(pH)x(total concentration)x(temperature)x(compounds)_, (if _rAllConc = True_).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Where species: [B], [B<sup>-</sup>], [B<sup>-2</sup>], [B<sup>-3</sup>]. 

### ThEq.solubilityHenry &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ThEq.solubilityHenry(compounds, wType='FW', temperature=[])
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature : _float_ or _ndarray of floats_, _optional, default: empty_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperature values.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Hs : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Henry's law solubility constant(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Shape: _(temperature)x(compounds)_, (if _temperature != empty_).<br>
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
ThP.getKeq(compounds, mRxn, temperature, phase)
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of the equilibrium constants.<br>

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

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

#### <ins>Bio-Thermodynamic State Analysis (BioThSA)</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Ecosystem modelling</ins>
:construction: Coming soon...

**·** _<ins>Kinetics</ins>_

The main and auxiliary functions to calculate kinetic rates are located in `reactions.py` module, and organized in two classes:
- `KinP`. Class with the functions to get the kinetic parameters from local databases.
- `KinRates`. Class with the functions to calculate the kinetic rates (biotic and abiotic transformations).

### KinP.getKinP &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
KinP.getKinP(typeParam, params, reaction, sample='All', comp=None)
```
Load the kinetic parameters from local databases (`kinetics\` folder).<br>
Return a dictionary with the requested kinetic parameters and an array with the `sample` names (rows of `typeParam.csv` file).<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeParam : _str_** <br>
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
KinRates.getRs(typeKin, typeParam, Rxn, Ct, sample='All', T=None, Tcorr=None)
```
Compute reaction rates as a function of limiting substrate, inhibitors and temperatures based on the chosen kinetic equation (`typeKin`) and temperature correlation (`Tcorr`).<br> 
Return a n-dimension array with the calculated kinetic rates and an array with the `sample` names (rows of `typeParam.csv` file).<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeKin : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Type of kinetic equations (e.g., rMM - 'Michaelis-Menten equation').<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **typeParam : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of parameter database, matching with csv name.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Rxn : _str_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested reaction.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Ct : _float, list or dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentration of substrates, products and/or inhibitors.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If Ct is a _dict_: {'Name compound': [concentrations]}.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **sample : _str or list_, _optional, default: 'All'_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Requested samples (rows of `typeParam.csv`).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **T : _float or list_, _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperatures.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Tcorr : _str_ ('Arrhenius'), _optional, default: None_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Resultant rates.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If Ct is a _dict_: {'Name compound': [concentrations]}.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rs : _np.ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Resultant substrate uptake rates.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If T and Tcoor are given: rs.shape = (parameters)x(concentrations)x(temperatures).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If they are not given: rs.shape = (parameters)x(concentrations).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **sampleNames : _str or list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of samples (_i.e._, rows of `typeParam.csv`).<br>

### KinRates.rMM &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
KinRates.rMM(qmax, Km, C)
```
Compute reaction rates as a function of limiting substrate  and inhibitors based on the Michaelis-Menten equation.<br> 
Return a n-dimension array with the calculated kinetic rates.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **qmax : _float, list, or np.ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of maximum uptake rate.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Km : _float, list, or np.ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Values of half-saturation constants (or Michaelis constants).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **C : _float, list, np.ndarray or dict_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Concentrations of limiting substrates, inhibitors or products.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If C is a _dict_: {'Name compound': [concentrations]}.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **r : _list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Resultant substrate uptake rates.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; r.shape = (parameters)x(concentrations).<br>

### KinRates.arrhCorr &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
KinRates.arrhCorr(rateBase, O, tempBase, temp)
```
Compute reaction rates as a function of temperature based on the Arrhenius correlation.<br> 
Return a n-dimension array with the calculated kinetic rates.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rateBase : _float, list, or np.ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Reaction rate at base (measured) temparature (tempBase).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **O : _float_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Arrhenius coefficient.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **tempBase : _float or list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; emperature base, that is, the original temperature of substrate rate.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temp : _float or list_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set of temperatures.<p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **rT : _np.ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Resultant substrate uptake rates.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; r.shape = (parameters)x(concentrations)x(temperatures).<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

## :clipboard: Instructions to use EcoSysEM platform via Command Line Interface (CLI)
:construction: Coming soon...

[🔼 Back to **Contents**](#readme-contents)
<!--
1. Download .zip code. Last version: `v#.#`. [Download package](https://github.com/soundslikealloy/EcoSysEM).
2. Extract files to a destination (Recommendation - Desktop).
3. Open **Anaconda Prompt or Terminal**.
4. Go to the **Code folder<sup>2</sup>** using `cd` command (more info about [Using Terminal](https://docs.anaconda.com/ae-notebooks/user-guide/basic-tasks/apps/use-terminal/?highlight=Using%20Terminal)).
    &#09;<br><sup><sup>2</sup>Code folder: folder with `ecosysem_cmd.py` file (Folder: `EcoSysEM\ecosysem`). </sup>
5. _Lorem ipsum..._
6. Execute one of the **EcoSysEM** blocks/functions using the following command lines.
-->

## Function Navigation
#### · <ins>Ideal Earth's atmosphere (ISA)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.setComposition](#isasetcomposition---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.selectRegion](#isaselectregion---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.getDictConc](#isagetdictconc---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.plotTandP](#isaplottandp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.plotCompsProfiles](#isaplotcompsprofiles---back-to-function-navigation)<br>

#### · <ins>MERRA2</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.getDataMERRA2](#merra2getdatamerra2---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.combDataMERRA2](#merra2combdatamerra2---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.selectRegion](#merra2selectregion---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.dictMERRA2](#merra2dictmerra2---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.keysMERRA2](#merra2keysmerra2---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [MERRA2.deleteKeyMERRA2](#merra2deletekeymerra2---back-to-function-navigation)<br>

#### · <ins>Thermodynamic equilibrium (ThEq)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThEq.plotpHSpeciation](#theqplotphspeciation---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThEq.pHSpeciation](#theqphspeciation---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThEq.solubilityHenry](#theqsolubilityhenry---back-to-function-navigation)<br>

#### · <ins>Thermodynamic parameters (ThP)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getThP](#thpgetthp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getDeltaG0r](#thpgetdeltag0r---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getDeltaH0r](#thpgetdeltah0r---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ThP.getKeq](#thpgetkeq---back-to-function-navigation)<br>

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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [KinRates.rMM](#kinratesrmm---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [KinRates.arrhCorr](#kinratesarrhcorr---back-to-function-navigation)<br>

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
