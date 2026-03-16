# SANDY - Installation and configuration

This guide explains how to install **SANDY** on **Windows** or **Linux** using pip, conda, or source code, and how to configure **NJOY** for full functionality (PENDF/ACE processing and random sampling).

## ⭐ 1. Installing SANDY with pip (Recommended)

### Prerequisites
 - Python installed (preferably via Anaconda).
 - Optional (but strongly recommended): NJOY2016 for nuclear data processing.

### Install
- Open a terminal and run:

```sh
# Upgrade pip to the latest version
python -m pip install --upgrade pip

# Install the latest stable version of SANDY from PyPI
pip install sandy
```

SANDY will install successfully even without NJOY, but some advanced features will be unavailable.

## ⭐ 2. Configuring NJOY (Optional but Recommended)

SANDY acts as a wrapper for the **NJOY** nuclear data processing code, enabling the generation of **PENDF** and **ACE** files.
While NJOY is not required for basic functionality, it is recommended, especially if you plan to produce random samples.

SANDY can automatically use the NJOY executable if you set the environment variable NJOY.

### LINUX
```sh
export NJOY=/path/to/njoy
```
To make this permanent, add it to your ``~/.bashrc``.

### Windows (Command Prompt)
```bat
set NJOY=C:\path\to\njoy.exe
```

### Verify NJOY detection
Open Python and run:
```python
import sandy
sandy.get_njoy()
```
If correctly configured, this prints the path to the NJOY executable.

## ⭐ 3. Creating a dedicated ``conda`` environment for sandy (Recommended)

Using a dedicated ``conda`` environment is the recommended approach to avoid dependency conflicts and keep your Python setup clean.
Both [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) and [Anaconda](https://www.anaconda.com/) provide a robust package manager and ship with many essential scientific Python libraries.

- **Miniconda** → a minimal installer containing only conda and Python
- **Anaconda** → a larger distribution that includes many pre‑installed scientific packages

### Create a clean environment (``sandy-devel``)
Once Miniconda/Anaconda is installed, open a terminal and create a new environment called ``sandy-devel``
```sh
conda update --name base conda
conda create -y --name sandy-devel -c conda-forge python numpy scipy pandas pyyaml pytables
```
Optional recommended packages:
 - **Data analysis**: `matplotlib`, `seaborn`, `scikit-learn`
 - **Testing**: `pytest`, `numpydoc` , `nbval`, `codecov`, `coveralls`, `pytest-cov`
 - **Packaging**: `build`, `twine`
 - **Notebooks**: `jupyterlab`, `jupyter_nbextensions_configurator`, `jupyter_contrib_nbextensions`

### Activate / deactivate environment
```sh
conda activate sandy-devel
conda deactivate
```
To manage your python environments read the [conda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).

### Install SANDY inside the conda environment
```sh
conda activate sandy-devel
pip install sandy
```

## ⭐ 4. Running SANDY in a python shell

1. Open a terminal or Anaconda Prompt.
2. Activate the environment (if using conda):
```sh
conda activate sandy-devel
``` 
3. Launch python (or ``ipython`` if installed):
```sh
python
``` 
5. Import SANDY:
```python
import sandy
```

If NJOY is configured, SANDY is fully ready.


## ⭐ 5. Installing SANDY from source
The source code of SANDY is available as a git repository. The recommended way to get SANDY is by cloning the source package using your local git distribution (*click [here](https://git-scm.com/downloads) to download it!*).

### Clone the repository
```sh
git clone https://github.com/luca-fiorito-11/sandy.git
cd sandy
```
If you encounter a **SSL certificate error**, disable SSL verification:
```git
git config --global http.sslVerify false
```

### Install from source
Inside the cloned sandy folder:
```sh
pip install .
```
(Activate your virtual environment beforehand, if using one.)

SANDY is installed! Now go back to run it in a python shell.


## ⭐ 6. Running SANDY with NJOY
SANDY also works as a wrapper to the NJOY code to process nuclear data files into **PENDF** and **ACE** formats.
The installation of NJOY is not mandatory (it is if you want to produce random samples), but recommended.

We recommend **NJOY2016** (NJOY2021 not yest tested), which can be found here:

👉 [https://github.com/njoy/NJOY2016](https://github.com/njoy/NJOY2016).

### NJOY on Linux

1. Clone the NJOY2016 repository.
2. Follow build instructions from the [NJOY documentation](https://docs.njoy21.io/install.html)
3. Export the NJOY path:
```sh
export NJOY=/path/to/njoy
```
4. Verify inside python:
```sh
import sandy
sandy.get_njoy()
```

### NJOY on Windows
NJOY requires a Linux-like environment. We recommend **Cygwin64**.
#### 1. Install Cygwin64
- Download [Cygwin64](https://cygwin.com/install.html).
- Follow the instructions of the installation wizard.
- You will be asked to select a 'Root Install Directory', that is, the directory where you want to install cygwin. In my case it is `C:\cygwin64
`. From now on we'll call the 'Root Install Directory' `C:\path\to\cygwin64`.
- Make sure you select the following packages to ensure that NJOY be succesfully installed:
    * `cmake 3.20.0-1`
    * `make 4.3-1`
    * `gcc-fortran 10.2.0-1`
    * `gcc-g++ 10.2.0-1`

#### 2. Download NJOY2016
- From a **git** terminal:
```sh
cd C:\path\to\cygwin64\home\your_username
git clone https://github.com/njoy/NJOY2016.git
```

#### 3. Build NJOY in Cygwin
- Open a `cygwin64` terminal and install NJOY2016:
```sh
cd C:\path\to\cygwin64\home\username\NJOY2016
mkdir bin
cd bin
cmake ..
make
make test
```
>  Make sure cmake finds an available python3 interpreter, if not you might have to use the cmake option `-DPython3_EXECUTABLE`.

#### 4. Register NJOY inside the conda environment
- Open an Anaconda Prompt terminal and set up the NJOY executable in the environment variable `NJOY`. This way SANDY will automatically find it.
```dos
conda activate sandy-devel
conda env config vars set NJOY=C:\path\to\cygwin64\home\username\NJOY2016\bin\njoy.exe
conda activate sandy-devel
```
- Check:
```dos
conda env config vars list
```

#### 5. Add Cygwin DLLs to PATH
To succesfully run NJOY Windows must be able to find some DLL files such as `cygwin1.dll`.

This file is part of cygwin, so most likely it's located in `C:\path\to\cygwin64\bin`.

Then, you have to add `C:\path\to\cygwin64\bin` (or the location where `cygwin1.dll` can be found) to your `PATH` typing the following on an Anaconda Prompt terminal:
```dos
set PATH=%PATH%;C:\path\to\cygwin64\bin
```

#### 6. Allow Cygwin to access Windows drives
If you want to succesfully run NJOY2016 through SANDY, cygwin must be allowed to access different directories outside the cygwin home directory.

From a cygwin terminal this can be done specifing `/cygdrive/` before the directory absolute path.

For example, you can access the root of your C: drive from cygwin by specifying the ``/cygdrive`` prefix:
```sh
cd /cygdrive/c
```
For convenience — for instance, to consistently be able to write in your user account on the C: drive (`C:\Users\your_username`) — create a symbolic link:
```sh
ln -sv /cygdrive/c/Users/your_username ~/your_username
```

#### 7. Test NJOY
From an Anaconda Prompt:
```sh
C:\path\to\cygwin64\home\your_username\NJOY2016\bin\njoy.exe
```
If NJOY runs, SANDY will be able to use it.

#### 8. Test NJOY in SANDY
Run python from an Anaconda Prompt:
```python
import sandy
sandy.get_njoy()
```

## ⭐ 7. Testing SANDY (only for source installations)
- Install required test packages:
```sh
conda install -y --name sandy-devel -c conda-forge pytest numpydoc nbval
```

- Run tests from an Anaconda Prompt from the ``sandy`` folder with:
```sh
conda activate sandy-devel
mkdir tests
cd tests
pytest ../sandy
```
All tests should pass successfully (this could take some time).

## ⭐ 8. Using SANDY in Jupyter Notebooks
For combatibility issues we recommend installing a python kernel specific for the `sandy-devel` environment.
For that, you can run the following after making sure that `ipykernel` is installed in the virtual environment.

```sh
conda activate sandy-devel
python -m ipykernel install --user --name sandy-devel --display-name "Python3 (sandy-devel)"
```

You can now select this kernel inside Jupyter.
