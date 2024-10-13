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
    - How to modify reaction databases | [GO](#how-to-modify-reaction-databases)
- Instructions to use EcoSysEM platform via Spyder | [GO](#how-to-modify-reaction-databases)
- EcoSysEM user guide | [GO](#ecosysem-user-guide)
    - EcoSysEM package layout | [GO](#ecosysem-package-layout)
    - Fundamentals and usage | [GO](#fundamentals-and-usage)
        - Environment creation and instance calling | [GO](#environment-creation-and-instance-calling)
        - Thermodynamic State Analysis (ThSA) | [GO](#thermodynamic-state-analysis-thsa)
        - Bio-Thermodynamic State Analysis (BioThSA) | [GO](#bio-thermodynamic-state-analysis-biothsa)
        - Ecosystem modelling | [GO](#ecosystem-modelling)
-  Instructions to use EcoSysEM platform via Command Line Interface (CLI) | [GO](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli)
-  Contact | [GO](#contact)
____________________________

## Before having fun...
**:warning: To open the links in a new tab: right click on the link + "Open link in new tab".**

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
[🔼 Back to **Contents**](#readme-contents)

____________________________

## :clipboard: Instructions for downloading and setting up EcoSysEM platform
1. Download .zip code. Last version: `v#.#`.<!-- [Download package](https://github.com/soundslikealloy/EcoSysEM). -->
2. Extract files to a destination (:bulb: Recommendation - Desktop).
3. Modify (if necessary) parameter databases using Excel files in folder  `\ecosysem\db\Excels` or built-in functions (see [How to modify parameter databases](#how-to-modify-parameter-databases) section).
4. Modify existing reaction databases or create a new one using Excel files in folder `\ecosysem\reactions\Excels` or built-in functions (see [How to modify reaction databases](#how-to-modify-reaction-databases) section).
5. Execute **EcoSysEM platform** via Spyder (see [Instructions to use EcoSysEM platform via Spyder](#clipboard-instructions-to-use-ecosysem-platform-via-spyder) section) or Command Line Interface (see [Instructions to use EcoSysEM platform via Command Line Interface (CLI)](#clipboard-instructions-to-use-ecosysem-platform-via-command-line-interface-cli) section).

[🔼 Back to **Contents**](#readme-contents)

### How to modify parameter databases
_Lorem ipsum..._

[🔼 Back to **Contents**](#readme-contents)

### How to modify reaction databases
_Lorem ipsum..._

[🔼 Back to **Contents**](#readme-contents)

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

[🔼 Back to **Contents**](#readme-contents)

### EcoSysEM package layout
Important modules and how to import functions or classes from them are listed below. Classes names start with a capital letter, functions with a lower letter, and attributes with a dot (.) and lower letter:
```
from ecosysem.module import function
from ecosysem.module import Class

EcoSysEM
  ├── envdef
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
  │           ├── .H2O
  │           ├── getVerticalProfiles
  │           ├── getDictConc
  │           ├── plotTandP
  │           └── plotCompsProfiles
  ├── thermodynamics
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
  │           └── plotDeltaGr
  ├── reactions
  │      └── Reactions
  │            ├── getRxn
  │            ├── getRxnByComp
  │            └── getRxnByName
  ├── ecosysem_spyder.py (Run EcoSysEM using Spyder, i.e., coding)
  └── ecosysem_cmd.py (Run EcoSysEM using Command Line Interface, CLI)
```

[🔼 Back to **Contents**](#readme-contents)

### Fundamentals and usage
_Lorem ipsum..._

[🔼 Back to **Contents**](#readme-contents)

#### <ins>Environment creation and instance calling</ins>
_Lorem ipsum..._

[🔼 Back to **Contents**](#readme-contents)

#### <ins>Thermodynamic State Analysis (ThSA)</ins>
_Lorem ipsum..._

[🔼 Back to **Contents**](#readme-contents)

#### <ins>Bio-Thermodynamic State Analysis (BioThSA)</ins>
:construction: Coming soon...

[🔼 Back to **Contents**](#readme-contents)

#### <ins>Ecosystem modelling</ins>
:construction: Coming soon...

[🔼 Back to **Contents**](#readme-contents)

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

[🔼 Back to top](#ecosysem-platform) &nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;[🔼 Back to **Contents**](#readme-contents)
