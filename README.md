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
- **NumPy**. NumPy is the fundamental package for scientific computing in Python. It is a Python library that provides a multidimensional array object, various derived objects (such as masked arrays and matrices), and an assortment of routines for fast operations on arrays, including mathematical, logical, shape manipulation, sorting, selecting, I/O, discrete Fourier transforms, basic linear algebra, basic statistical operations, random simulation and much more. For more info and tutorials, click [here](https://numpy.org/).
- **Pandas**. Pandas is a fast, powerful, flexible and easy to use open source data analysis and manipulation tool, built on top of the Python programming language. For more info and tutorials, click [here](https://pandas.pydata.org/).
- **Matplotlib**. Matplotlib is a library for creatinc static, animated and interactive visualizations in Python. For more info and tutorials, click [here](https://matplotlib.org/).
- **SciPy**. SciPy is a collection of mathematical algorithms and convenience functions built on NumPy . It adds significant power to Python by providing the user with high-level commands and classes for manipulating and visualizing data. For more info and tutorials, click [here](https://scipy.github.io/devdocs/tutorial/index.html).
- **iteration_utilites**. Iteration_utilities is a collection of functional programming based on and utilizing iteratiors and generators. Most of the functions are inspiered by the _itertools_ module, but implemented in C to achieve a better overall performance. For more info and tutorials, click [here](https://iteration-utilities.readthedocs.io/).

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
**iteration_utilites**:
```
pip install iteration_utilities
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

|  Compound  |  Reaction no 1.0  |  Reaction no 2.0  |  ...  |  Reaction no N.0  |
| ---------- | ----------------- | ----------------- | ----- | ----------------- |
| Compound 1 | Stoich. value 1.1 | Stoich. value 2.1 |  ...  | Stoich. value N.1 |
| Compound 2 | Stoich. value 1.2 | Stoich. value 2.2 |  ...  | Stoich. value N.1 |

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
  │      ├── Environment
  │      │    ├── .temperature
  │      │    ├── .pressure
  │      │    ├── .pH
  │      │    ├── .compounds
  │      │    ├── .compositions
  │      │    ├── setT
  │      │    ├── setP
  │      │    ├── setpH
  │      │    └── setComposition
  │      └── ISA {subclass of Environment}
  │           ├── ._ISAproperties
  │           ├── .dryComposition
  │           ├── .altitude
  │           ├── .temperature
  │           ├── .pressure
  │           ├── .compounds
  │           ├── .compositions
  │           ├── getVerticalProfiles
  │           ├── getDictConc
  │           ├── plotTandP
  │           └── plotCompsProfiles
  ├── reactions.py 
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
  - General Environment | [GO](#general-environment)
  - Ideal Earth's atmosphere (International Standard Atmosphere, ISA) | [GO](#ISA)
  <!-- - How to create a new Environment subclass | [GO](#create-new-environment-subclass) -->
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

<a name="general-environment">**General environment**</a><br>
Environments can be defined as class instances. An instance is an object that's built from a class and contains real data. Many instances can be created from a single class. To create a new general environment instance (_i.e.,_ instantiate the class `Environment`), several instances attributes must be given, called `.temperature`, `.pressure`, `.pH`, `.compounds` and `.compositions`:
```python
from envdef import Environment

newEnv1 = Environment(-5, 5.0, 7.0, ['A', 'B', 'C'], [1.00e-3, 1.00e-3, 1.00e-3])
newEnv2 = Environment(15, 1.0, 7.0, ['O2', 'CO2', 'CH4', 'NH3'], [[9.375e-5, 1.0e-7], [1.000e-3, 1.000e-3], [3.500e-3, 2.231e-5], [1.000e-3, 1.000e-5]])
```
Once created the Environment instances (`newEnv1` and `newEnv2`), their attributes can be accessed using **dot notation**. All attributes of `Environment` class are shown in [EcoSysEM package layout](#ecosysem-package-layout). Here are some examples:
```python
>>> print(newEnv1.temperature)
-5
>>> print(newEnv2.temperature)
15
>>> print(newEnv1.compositions)
{'A': 0.001, 'B': 0.001, 'C': 0.001}
>>> print(newEnv2.compositions)
{'O2': [9.375e-05, 1e-07], 'CO2': [0.001, 0.001], 'CH4': [0.0035, 2.231e-05], 'NH3': [0.001, 1e-05]}
```
In addition to class attributes, `Environment` class also contains class functions (known as _instance mehtods_). These functions can only be called on an instance of that class. All functions of `Environment` class are shown in [EcoSysEM package layout](#ecosysem-package-layout). The `Environment` class has four instance methods:

### Environment.setT &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.setT(newT)
```
Modify temerature of the `Environment` instance.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **newT : _int_ or _float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; New tempearture value of `Environment` instance.

### Environment.setP &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.setP(newP)
```
Modify pressure of the `Environment` instance.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **newP : _int_ or _float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; New pressure value of `Environment` instance.

### Environment.setpH &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.setpH(newpH)
```
Modify pH of the `Environment` instance.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **newpH : _int_ or _float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; New pH value of `Environment` instance.

### Environment.setComposition &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Environment.setComposition(compound, composition)
```
Add (if _compound_ does not exist) or modify (if _compound_ does exist) composition of the `Environment` instance.<p>
**Parameters:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **compound : _str_ or _list of strs_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; New compound or compound to modify.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **composition : _float_ or _list of float_**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Composition of new or existing compound.<br>

To call the instance methods of `Environment`, the name of the `Environment` object must be preceded by the function name and its arguments in parentheses. Here is an example:
```python
from envdef import Environment

newEnv = Environment(-5, 5.0, 7.0, ['A', 'B', 'C'], [1.00e-3, 1.00e-3, 1.00e-3])

>>> print(newEnv.pressure)
5.0
>>> print(newEnv.pH)
7.0
>>> print(newEnv.compositions)
{'A': 0.001, 'B': 0.001, 'C': 0.001}

newEnv.setP(1.0)
newEnv.setpH(8.5)
newEnv.setComposition(['A', 'D'], [5.00e-3, 2.00e-3])

>>> print(newEnv.pressure)
1.0
>>> print(newEnv.pH)
8.5
>>> print(newEnv.compositions)
{'A': 0.005, 'B': 0.001, 'C': 0.001, 'D': 0.002}
```

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

From `Environment` class, new inhereted classes (also known as _subclasses_) can be created inhereting the attributes. New attributes and function can be defined in these sublcasses, which will belong only to the subclass in question.

<a name="ISA">**Ideal Earth's atmosphere (International Standard Atmosphere, ISA)**</a><br>
The International Standard Atmosphere (ISA) is a static atmospheric model of how pressure, temperature, density and viscosity of the Earth's atmosphere change over a a wide range of altitudes. The ISA model is detailed in ISO 2533:1975[^1]. The ISA model divides the atmosphere into different layers with specific physical properties (such as temperature rate change, base temperature, base atmospheric pressure or base atmospheric density). With these properties, the temperature and pressure profiles are calculated. On the other hand, the ISA model assumes that Earth's atmosphere is a mixture of gas, water vapor and a certain quantity of aerosols, in which the composition remains pratically constant up to altitudes of 90 - 95 km. The dry (0%<sub>vol</sub> of water) and wet composition (up to 4%<sub>vol</sub> of water) of Earth's atmosphere are from National Oceanic and Atmospheric Administration (NOAA)[^2] and _Atmospheric Radiation: Theoretical Basis (2<sup>nd</sup> edition)_[^3] (**Table 1**).

**Table 1. Dry Earth's atmosphere composition (in %<sub>vol</sub>)**

| Compound | Value | Compound | Value | Compound | Value |
|---|---|---|---|---|---|
| __N<sub>2</sub>__  | 7.8084·10<sup>-1</sup> | __H<sub>2</sub>__  | 5.500·10<sup>-7</sup> | __NO__             | 1.000·10<sup>-9</sup>  |
| __O<sub>2</sub>__  | 2.0946·10<sup>-1</sup> | __N<sub>2</sub>O__ | 3.300·10<sup>-7</sup> | __SO<sub>2</sub>__ | 1.000·10<sup>-9</sup>  |
| __Ar__             | 9.3400·10<sup>-3</sup> | __CO__             | 1.000·10<sup>-7</sup> | __H<sub>2</sub>S__ | 5.000·10<sup>-11</sup> |
| __CO<sub>2</sub>__ | 4.2000·10<sup>-4</sup> | __Xe__             | 9.000·10<sup>-8</sup> | | |
| __Ne__             | 1.8182·10<sup>-5</sup> | __O<sub>3</sub>__  | 7.000·10<sup>-8</sup> | | |
| __He__             | 5.2400·10<sup>-6</sup> | __NO<sub>2</sub>__ | 2.000·10<sup>-8</sup> | | |
| __CH<sub>4</sub>__ | 1.9200·10<sup>-6</sup> | __I<sub>2</sub>__  | 1.000·10<sup>-8</sup> | | |
| __Kr__             | 1.1400·10<sup>-6</sup> | __NH<sub>3</sub>__ | 4.000·10<sup>-9</sup> | | |

To create a new _ISA_ object (_i.e.,_ instantiate the subclass `ISA`), the instance attributes `layers`, `H2O`, `pH` and `resolution` are necessary:
- `layers`. Selection of atmosphere layers defined by ISA model[^1]. This attribute can be 'All' (_string_), an _integer_ from 0 to 7 or a _list_ of integers.
- `H2O`. Water content of atmosphere. This attribute must be a _float_ from 0.0 to 0.04.
- `pH`. pH of atmosphere. This attribute must be a _float_.
- `resolution`. Resolution of altitude array, that is, the size of altitude nodes per layer (in m). This attribute must be an _integer_.

Because `ISA` sublcass is a inhereted class of `Environment` class, this has also `.temperature`, `.pressure`, `.pH`, `.compounds` and `.composistions`. Additionally, `ISA` subclass has also its own inherent attributes, that is, attributes that are part of the essential nature of `ISA` sublcass: _i)_ the properties of ISA layers (`._ISAproperties`), _ii)_ dry composition (`.dryComposition`), _iii)_ altitude (an NumPy array with the range of altitudes of ISA instance).
`ISA` subclass also contains its own class functions (or _instance mehtods_). All functions of `ISA` subclass are summarized in [EcoSysEM package layout](#ecosysem-package-layout). The `ISA` subclass has four instance methods:

### ISA.getVerticalProfiles &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.getVerticalProfiles(phase, compound=None)
```
Return vertical profiles in _format=ndarray_ of selected compounds or all compounds of `ISA` subclass.<p>
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Pi, concG, compoundsG : _ndarray_** (if _phase='G'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **concLFW, compoundsLFW : _ndarray_** (if _phase='L-FW'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **concLSW, compoundsLSW : _ndarray_** (if _phase='L-SW'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **concLFW, concLSW, compoundsLFW, compoundsLSW : _ndarray_** (if _phase='L'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Pi, concG, concLFW, concLSW, compoundsG, compoundsLFW, compoundsLSW  : _ndarray_** (if _phase='All'_)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Pi : Parcial pressure of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; concG : Concentration in gas of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; concLFW : Concentration in liquid (freshwater) of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; concSFW : Concentration in liquid (seawater) of desired compounds.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; compoundsG : Compounds names (gas).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; compoundsLFW : Compounds names (freshwater).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; compoundsLFW : Compounds names (seawater).<br>

### ISA.getDictConc &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.getDictConc(phase, compound=None)
```
Return vertical profiles in _format=dict_ of selected compounds or all compounds of `ISA` subclass.<p>
**Parameters:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **phase : _str ('G', 'L-FW', 'L-SW', 'L' or 'All')_ **<br>
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

### ISA.plotTandP &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.plotTandP()
```
Plot temperature and pressure profiles of `ISA` instance.<p>
**Parameters:** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **None** <p>
**Returns:**<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Spyder plot** <br>

### ISA.plotCompsProfiles &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
ISA.plotCompsProfiles(Conc, xLabel, logCLabel=False, compounds=None)
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

<!--
<a name="create-new-environment-subclass">**How to create a new Environment subclass**</a><br>
_Lorem ipsum_

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

-->

#

#### <ins>Ecosystem Analysis (EcoA)</ins>
:construction: Coming soon...

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

#### <ins>Thermodynamic State Analysis (ThSA)</ins>
With the <ins>Thermodynamic State Analysis module (ThSA)</ins>, the user can identify which rections are feasible (_i.e.,_ exergonic) considering the environmental conditions (nonstandard conditions). It is important to take into account that because a given reaction is exergonic under a particular environmental conditions (ΔG<sub>r</sub><0), it does not necessary mean that organisms will be able to catalyze it and, if they can, they have enough energy for growth and/or maintenance. Remember that Gibbs energy quantifies the tendency of a chemical reaction to proceed in a particular direction. To determine the actual energy that organisms can take from a particular transformation and if these can growth and/or satisfy maintenance requirements, the [Biological Thermodynamic State Analysis (BioThSA)](#bio-thermodynamic-state-analysis-biothsa) should be used. To calculate the nonstandard Gibbs free energy of reaction (ΔG<sub>r</sub>), the Gibbs-Helmholtz-Nernst relationship is used. The Gibbs-Helmholtz-Nernst relationship considers the influnce of temperature, pH and concentrations of substrates and products on ΔG<sub>r</sub>.

The main and auxiliary functions to perform the ThSA are located in `thermodynamics.py` module, and organized in three classes:
- `ThSA`. Class with the functions to perform the ThSa and export the results.
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **temperature : _float_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Temperature value.<p>
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user. The function always returns the abbreviation.<br>

### Reactions.getRxnByComp() &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user. The function always returns the abbreviation.<br>

### Reactions.getRxnByName() &nbsp;&nbsp;&nbsp;&nbsp; <sup><sub>[🔽 Back to Function Navigation](#function-navigation)</sub></sup>
```python
Reactions.getRxnByName(typeRxn, nameRxn, warnings = False)
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Stoichiometric matrix of reactions.<p>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **infoRxn : _ndarray_** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name of reactions given by the user. The function always returns the abbreviation.<br>

[🔼 Back to **Fundamentals and usage**](#fundamentals-and-usage) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

#

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
2. Extract files to a destination (Recommendation - Desktop).
3. Open **Anaconda Prompt or Terminal**.
4. Go to the **Code folder<sup>2</sup>** using `cd` command (more info about [Using Terminal](https://docs.anaconda.com/ae-notebooks/user-guide/basic-tasks/apps/use-terminal/?highlight=Using%20Terminal)).
    &#09;<br><sup><sup>2</sup>Code folder: folder with `ecosysem_cmd.py` file (Folder: `EcoSysEM\ecosysem`). </sup>
5. _Lorem ipsum..._
6. Execute one of the **EcoSysEM** blocks/functions using the following command lines.
-->

## Function Navigation
#### · <ins>General environment (_for any subclass of Environment_)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.setT](#environmentsett---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.setP](#environmentsetp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.setpH](#environmentsetph---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Environment.setComposition](#environmentsetcomposition---back-to-function-navigation)<br>

#### · <ins>Ideal Earth's atmosphere (ISA)</ins>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.getVerticalProfiles](#isagetverticalprofiles---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.getDictConc](#isagetdictconc---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.plotTandP](#isaplottandp---back-to-function-navigation)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [ISA.plotCompsProfiles](#isaplotcompsprofiles---back-to-function-navigation)<br>

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

[🔼 Back to **Contents**](#readme-contents)
__________________________________________________

## Contact

**Eloi Martinez-Rabert**. :envelope: eloi.mrp@gmail.com

[🔼 Back to **Top**](#ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)

### References

[^1]: International Organization for Standardization, Standard Atmosphere, ISO 2533:1975, 1975. [Link to ISO](https://www.iso.org/standard/7472.html).
[^2]: The Atmosphere. Introduction to the Atmosphere. From National Oceanic and Atmospheric Administration (NOAA). Last updated: 2 July, 2024 [Link to NOAA](https://www.noaa.gov/jetstream/atmosphere).
[^3]: Goody, R. M., & Yung, Y. L. (1989). Atmospheric Radiation: Theoretical Basis. Oxford University Press. [https://doi.org/10.1093/oso/9780195051346.001.0001](https://doi.org/10.1093/oso/9780195051346.001.0001).
