# Installation

## Download

The source code of SANDY is available as a git repository. The recommended way to get SANDY is by cloning the source package using your local git distribution (*click [here](https://git-scm.com/downloads) to download it!*).

```sh
git clone https://github.com/luca-fiorito-11/sandy.git
cd sandy
```

If a `SSL Certificate problem` appears when pushing or pulling with git, one can tell git to not perform the validation of the certificate using the global option:

```git
git config --global http.sslVerify false
```

## Installation using Anaconda

### Downloading Anaconda
The recommended way to install the SANDY requirements (especially on Windows) is with the Anaconda distribution, which already includes several python packages and few other cool tools. The installer for the most recent Anaconda distributions can be found [here](https://www.anaconda.com/products/individual). For the Anaconda installation documentation, click [here](https://docs.anaconda.com/anaconda/install/).

### Getting started with Anaconda
The user guide to get yourself acquainted with Anaconda is available [here](https://docs.anaconda.com/anaconda/user-guide/). On a Windows machine we recommend that you use the Anaconda Prompt terminal.

### Installing SANDY

#### Creating a clean virtual environment using conda (not mandatory, but recommended)
We advise to install the SANDY dependences in a python environment that was not previously altered. You can do so from a terminal/Anaconda Prompt using the package/environment management system `conda` included in Anaconda, as
```sh
conda update --name base conda
conda create -y --name sandy-devel -c conda-forge python=3.7 jupyterlab jupyter_nbextensions_configurator
```
I like to include jupyterlab to my environment!

Then, dependences are installed from the requirement file,
```sh
conda install --force-reinstall -y --name sandy-devel -c conda-forge --file requirements_conda.txt
```
Make sure you are in the folder `sandy` that you cloned with git and where file `requirements_conda.txt` is located.

From now on, every time you want to use SANDY in your python distribution you need to activate the correct environemt.
```sh
conda activate sandy-devel
```
... and the following to deactivate it.
```sh
conda deactivate
```
To manage your python environments read the [conda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).

#### Working without a virtual environment
Move to the folder `sandy` that you cloned with git and install the requirements

**with conda**
```sh
conda install --force-reinstall -y -c conda-forge --file requirements_conda.txt
```

**with pip**
```sh
pip install -r requirements.txt
```

> The difference between files `'requirements.txt'` and `'requirements_conda.txt'` is all in how they call the python package `pytables`. For more info have a look at the pytables [installation website](https://www.pytables.org/usersguide/installation.html). If you are aware of a better way of doing it, please get in touch with us!

### Installation using `setup.py`
SANDY is still using some fortran (hopefully not for long!), which is compiled with the `f2py` module of `numpy`. For that an updated C++ compiler and a Fortran compiler are needed.

**on Linux**
We just assume the compilers mentioned above are available. Then, once you are in the folder `sandy`, run the installation command,
```sh
python setup.py install
```
> Don't forget to activate the virtual environment if you are using one.

**on Windows**

For the C++ compiler, [Microsoft C++ Build Tools]((https://visualstudio.microsoft.com/visual-cpp-build-tools/)) are recommended. Once the download has completed follow these steps:
1. Launch the Visual Studio Build Tools installation wizard.
2. In the `Workloads` tab select the `C++ build tools package` and click on `Install`.
3. A reboot of your machine is requested.

For the Fortran compiler do the following:
1. Open the `Anaconda Navigator` from the start menu.
2. Select `sandy-devel` or `base (root)` in the `Environments` tab, respectively if you created or not a virtual environment.
3. On the top widget change the selection from `Installed` to `All`.
4. In the `Search Packages` box type *fortran*.
3. Make sure that `m2w64-gcc-fortran` and `m2w64-gcc-libgfortran` are selected and then click on `Apply` on the bottom right-side.
4. Click once again on the `Apply` button in the pop-up window to accept package modifications (dependencies).

Once the C++ compiler and the Fortran compiler are installed, open an Anaconda Prompt terminal, activate the virtual environment (if any) and install SANDY as
```sh
conda activate sandy-devel
python setup.py install
```

Go to the `site-packages` directory containing the installed python packages. You can find this directory running the following inside a python shell:
```python
from distutils.sysconfig import get_python_lib
print(get_python_lib())
```
It should be something like `C:\path\to\Anaconda\envs\sandy-devel\Lib\site-packages` (or `C:\path\to\Anaconda\Lib\site-packages`if you don't work with a virtual environment).

There you should find a *python egg* for SANDY, that is, a folder named `sandy-*.egg`. If you move to folder `.\rwf\.libs` inside the SANDY egg you'll see that some `librwfortra.*.gfortran-win_*.dll` files were created during the installation. Copy such `*.dll` files at the root level of the SANDY egg (`C:\path\to\Anaconda\envs\sandy-devel\Lib\site-packages\sandy-*.egg` or `C:\path\to\Anaconda\Lib\site-packages\sandy-*.egg`)

Now SANDY should work also on Windows! If it still doesn't, some dependencies might be missing (you might see an import module error). In that case, try the following:
1. Download a [software](https://visualstudio.microsoft.com/visual-cpp-build-tools/) to scan module dependencies.
2. Open `DependenciesGui.exe` and load the `*.pyd` file of the module SANDY is failing to import. Missing dependencies (DLL files) will appear with a `?` sign.
3. Search for the missing dependencies on your PC and copy them to the same directory where the `*.pyd` file of the failing module is located.
4. Good luck!


### Adding SANDY to PYTHONPATH
If you want to add the `sandy` folder to your `PYTHONPATH` you can do

**in Linux**
```sh
export PYTHONPATH="/path/to/sandy:${PYTHONPATH}"
```
Add the previous command to your `${HOME}/.bashrc` or `${HOME}/.bash_profile` file to make your changes permanent.

**in Windows**
```dos
set PYTHONPATH=%PYTHONPATH%;C:\path\to\sandy
```
or, if you want the changes to be permanent,
```dos
setx PYTHONPATH=%PYTHONPATH%;C:\path\to\sandy
```

## Running a python shell
1. Open a terminal/Anaconda Prompt.
2. Switch to the correct python environment (if any).
3. Type `python` or `ipython` to open up a python/ipython session.
4. Then type
```python
import sandy
```

Now, if `sandy` is available you are ready to roll!

## Testing SANDY
A number (far away from being extensive enough) of unit tests is available to test that SANDY works correctly. First, make sure you istalled the following packages

**with Anaconda**
```sh
conda install -y --name sandy-devel -c conda-forge numpydoc nbval
```

**with pip**
```sh
pip install numpydoc nbval
```

Then, you can run the tests from a terminal/Anaconda Prompt from the `sandy` folder with
```sh
conda activate sandy-devel
pytest sandy
```

Hopefully they are all green!


## Running SANDY in a Jupyter notebook

For combatibility issues we recommend installing a python kernel specific for the `sandy-devel` environment.
For that, you can run the following after making sure that `ipykernel` is installed in the virtual environment.

```sh
conda activate sandy-devel
python -m ipykernel install --user --name sandy-devel --display-name "Python3 (sandy-devel)"
```

## Running SANDY with NJOY
SANDY also works as a wrapper to the NJOY code to process nuclear data files into PENDF and ACE formats.
The installation of NJOY is not mandatory, but recommended.

We suggest using NJOY2016 (we haven't tried NJOY2021 yet!).
You can find the source on its [github repository](https://github.com/njoy/NJOY2016).

To install NJOY you can do the following

**in Linux**

Clone NJOY2016 with `git` and follow the installation instructions provided on the [NJOY website](https://docs.njoy21.io/install.html).

**in Windows**

1. Download [Cygwin64](https://cygwin.com/install.html).
2. Follow the instructions of the installation wizard.
3. You will be asked to select a 'Root Install Directory', that is, the directory where you want to install cygwin. In my case it is `C:\cygwin64
`. From now on we'll call the 'Root Install Directory' `C:\path\to\cygwin64`.
3. Make sure you select the following packages to ensure that NJOY be succesfully installed:
    * `cmake 3.20.0-1`
    * `make 4.3-1`
    * `gcc-fortran 10.2.0-1`
    * `gcc-g++ 10.2.0-1`
4. Open a `git` terminal and download NJOY2016.
```sh
cd C:\path\to\cygwin64\home\your_username
git clone https://github.com/njoy/NJOY2016.git
```
5. Open a `cygwin64` terminal and install NJOY2016:
```sh
cd C:\path\to\cygwin64\home\username\NJOY2016
mkdir bin
cd bin
cmake ..
make
make test
```
>  Make sure cmake finds an available python3 interpreter, if not you might have to use the cmake option `-DPython3_EXECUTABLE`.
6. Open an Anaconda Prompt terminal and set up the NJOY executable in the environment variable `NJOY`. This way SANDY will automatically find it.
```dos
conda activate sandy-devel
conda env config vars set NJOY=C:\path\to\cygwin64\home\username\NJOY2016\bin\njoy.exe
conda activate sandy-devel
```

If you run the following you should see `NJOY` in your list of environment variables,
```dos
conda env config vars list
```

To verify that NJOY is correctly found by SANDY open a python terminal, import sandy and run `get_njoy()`.
```python
import sandy
sandy.get_njoy()
```

## `git` branches
The development of SANDY is mostly carried out on the `git` branch `develop`. In this branch you should have the latest features available. We recommend using this branch if you want to use SANDY as an interface to nuclear data files. All SANDY tutorials available [here](https://luca-fiorito-11.github.io/sandy_notebooks/) were produced with the `develop` branch and might not work on other branches. 

The `master` branch was originally developed for creating random perturbed files of evaluated nuclear data libraries. If this is your objective, that is, you could't care less of how SANDY works and just want some results, then you should use the `git` branch `master`.

To have a list of the available branches open a `git` shell, move to the `sandy` folder and type,
```git
git branch --list
```

If you want to switch to, say, the `master` branch type 
```git
git checkout master
```

If SANDY is already installed on your PC and you know that changes were made on the corresponding remote version on github you don't to clone the code once again. Just move to the `sandy` folder and do
```git
git pull
```
from the selected branch.

If any conflict appear, then it is probably time to read the `git` [documentation](https://git-scm.com/doc).
