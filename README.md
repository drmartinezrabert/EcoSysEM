# EcoSysEM platform
Eco-System Evolution &amp; Modelling

*Contributors: Eloi Martinez-Rabert*
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
        - Bio-Thermodynamic State Analysis (BioThSA) | [GO](#bio-thermodynamic-state-analysis-biothsa)
        - Ecosystem modelling | [GO](#ecosystem-modelling)
-  Instructions to use EcoSysEM platform via Command Line Interface (CLI) | [GO](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli)
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
- **NumPy**. NumPy is the fundamental package for scientific computing in Python. It is a Python library that provides a multidimensional array object, various derived objects (such as masked arrays and matrices), and an assortment of routines for fast operations on arrays, including mathematical, logical, shape manipulation, sorting, selecting, I/O, discrete Fourier transforms, basic linear algebra, basic statistical operations, random simulation and much more. For more info and tutorials, click [here](https://numpy.org/).
- **Pandas**. Pandas is a fast, powerful, flexible and easy to use open source data analysis and manipulation tool, built on top of the Python programming language. For more info and tutorials, click [here](https://pandas.pydata.org/).
- **Matplotlib**. Matplotlib is a library for creatinc static, animated and interactive visualizations in Python. For more info and tutorials, click [here](https://matplotlib.org/).
- **SciPy**. SciPy is a collection of mathematical algorithms and convenience functions built on NumPy . It adds significant power to Python by providing the user with high-level commands and classes for manipulating and visualizing data. For more info and tutorials, click [here](https://scipy.github.io/devdocs/tutorial/index.html).

[🔼 Back to **Contents**](#readme-contents)

### Installation of packages using Anaconda Navigator
You can install any Python package using the **Anaconda Navigator**. For this, execute the navigator and click to **Environments**. In this section you can install new packages and delete the already installed. For more info, click [here](https://docs.anaconda.com/free/navigator/).

[🔼 Back to **Contents**](#readme-contents)

### Installation of packages using pip
**pip** is the package installer for Python. In general, pip installs the minimal instalation requirements automatically, but not the optionals requirements. To install the mentioned packages using pip, you have only to write the following command lines in **Anaconda Prompt or Terminal**:

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
[🔼 Back to **Contents**](#readme-contents)

____________________________

## :clipboard: Instructions for downloading and setting up EcoSysEM platform
1. Download .zip code. Last version: `v#.#`.<!-- [Download package](https://github.com/soundslikealloy/EcoSysEM). -->
2. Extract files to a destination (:bulb: Recommendation - Desktop).
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
- metabolisms. Stoichiometrix matrix of metabolic reactions.

To modify existing databases open de .xlsx file in `ecosysem\reactions\Excels\*.xlsx`. The user can also create a new database reaction. Once finished, **1)** _Save_ Excel file in `\ecosysem\reactions\Excels\` folder and **2)** _Save as_ the document in .csv format in `\ecosysem\reactions\` folder. All reaction databases have the same structure:

|  Compound  |  Reaction no 1.0  |  Reaction no 2.0  |  ...  |  Reaction no N.0  |
| ---------- | ----------------- | ----------------- | ----- | ----------------- |
| Compound 1 | Stoich. value 1.1 | Stoich. value 2.1 |  ...  | Stoich. value N.1 |
| Compound 2 | Stoich. value 1.2 | Stoich. value 2.2 |  ...  | Stoich. value N.1 |

Where **Stoich. value A.B** is the stoichiometric value of _Compound B_ for _Reaction A_. Stoichiometric values are negative for substrates (< 0) and positive for products (> 0). If a compound does not participate in a reaction, that cell is left blank. Each column is a specific reaction, and in the headers In column headers, it is written the participating compounds between '/' and the name of reaction between parenthesis.  

· For example, biotic ammonia oxidation to nitrate by comammox bacteria (CMX): NH<sub>3</sub> + 2.0·O<sub>2</sub> → NO<sub>3</sub><sup>-</sup> + H<sub>2</sub>O + H<sup>+</sup>.

|  Compound  |  NH3/O2/NO3-/H2O/H+ (CMX)  |
| ---------- | -------------------------- |
| H+ | 1 |
| H2O | 1 |
| O2 | -2 |
| NH3 | -1 |
| NO2- |  |
| NO3- | 1 |

> [!NOTE]
> In pHSpeciation database, a specific reaction names (in parenthesis) are used:
> - First deP: first deprotonation.
> - Second deP: second deprotonation.
> - Third deP: third deprotonation.

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
This guide is an overview and explains the important features of **EcoSysEM platform**.

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### EcoSysEM package layout
Important modules and how to import functions or classes from them are listed below. Classes names start with a capital letter, functions with a lower letter, and attributes with a dot (.) and lower letter:
```python
from ecosysem.module import function
from ecosysem.module import Class

ecocysem
  ├── envdef.py 
  │      ├── Environment
  │      │      ├── .temperature
  │      │      ├── .pressure
  │      │      ├── .pH
  │      │      ├── .compounds
  │      │      ├── .compositions
  │      │      ├── setT
  │      │      ├── setP
  │      │      ├── setpH
  │      │      └── setComposition
  │      └── ISA {subclass of Environment}
  │           ├── .altitude
  │           ├── .temperature
  │           ├── .pressure
  │           ├── .compounds
  │           ├── .compositions
  │           ├── getVerticalProfiles
  │           ├── getDictConc
  │           ├── plotTandP
  │           └── plotCompsProfiles
  ├── thermodynamics.py 
  │      ├── ThP
  │      │    ├── getThP
  │      │    ├── getDeltaG0r
  │      │    ├── getDeltaH0r
  │      │    └── getKeq
  │      ├── ThEq
  │      │     ├── solubilityHenry
  │      │     ├── pHSpeciation
  │      │     └── plotpHSpeciation
  │      └── ThSA
  │           ├── getDeltaGr
  │           └── exportDeltaGr
  ├── reactions.py 
  │      └── Reactions
  │            ├── getRxn
  │            ├── getRxnByComp
  │            └── getRxnByName
  ├── ecosysem_spyder.py
  └── ecosysem_cmd.py
```

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### Fundamentals and usage
This section clarifies concepts, design decisions and technical constraints in **EcoSysEM platform**. **EcoSystem platform** constitutes on four main blocks:
- Environment definition | [GO](#environment-definition-and-instance-calling)
- Thermodynamic State Analysis (ThSA) | [GO](#thermodynamic-state-analysis-thsa)
- Bio-Thermodynamic State Analysis (BioThSA) | [GO](#bio-thermodynamic-state-analysis-biothsa)
- Ecosystem modelling | [GO](#ecosystem-modelling)

[🔼 Back to **Instructions (EcoSysEM via Spyder)**](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Environment definition and instance calling</ins>
_Lorem ipsum..._

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Thermodynamic State Analysis (ThSA)</ins>
_Lorem ipsum..._

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Bio-Thermodynamic State Analysis (BioThSA)</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#### <ins>Ecosystem modelling</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

## :clipboard: Instructions to use EcoSysEM platform via Command Line Interface (CLI)
:construction: Coming soon...

[🔼 Back to **Contents**](#readme-contents)
<!--
1. Download .zip code. Last version: `v#.#`. [Download package](https://github.com/soundslikealloy/EcoSysEM).
2. Extract files to a destination (:bulb: Recommendation - Desktop).
3. Open **Anaconda Prompt or Terminal**.
4. Go to the **Code folder<sup>2</sup>** using `cd` command (more info about [Using Terminal](https://docs.anaconda.com/ae-notebooks/user-guide/basic-tasks/apps/use-terminal/?highlight=Using%20Terminal)).
    &#09;<br><sup><sup>2</sup>Code folder: folder with `ecosysem_cmd.py` file (Folder: `EcoSysEM\ecosysem`). </sup>
5. _Lorem ipsum..._
6. Execute one of the **EcoSysEM** blocks/functions using the following command lines.
-->
____________________________

## Contact

**Eloi Martinez-Rabert**. :envelope: eloi.mrp@gmail.com

[🔼 Back to **Top**](#ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)
