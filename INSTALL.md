# Installation and configuration

## Installation on Windows/Linux with pip/conda
The recommended way to install SANDY both on Windows and Linux is on an Anaconda virtual environment that already includes several python packages. Click [here](#anaconda) for more details on Anaconda.

We advise to install the SANDY dependences in a python environment that was not previously altered. You can do so from a terminal/Anaconda Prompt using the package/environment management system `conda` included in Anaconda, as
```sh
conda update --name base conda
conda create -y --name sandy-devel -c conda-forge python numpy scipy pandas pyyaml pytables
```
This environment covers all hard dependencies of SANDY, but i usually like to add a few more packages (such as `ipython`) depending on my goal:
 - data analysis and visualization: `matplotlib`, `seaborn`, `scikit-learn`
 - testing: `pytest`, `numpydoc` , `nbval`, `codecov`, `coveralls`, `pytest-cov`
 - packaging projects: `build`, `twine`
 - jupyter notebooks: `jupyterlab`, `jupyter_nbextensions_configurator`, `jupyter_contrib_nbextensions`

From now on, every time you want to use SANDY in your python distribution you need to activate the correct environemt.
```sh
conda activate sandy-devel
```
... and the following to deactivate it.
```sh
conda deactivate
```
To manage your python environments read the [conda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).

Once the virtual environment is set up and active, you can install SANDY with `pip`.

```sh
pip install -i https://test.pypi.org/simple/ sandy
```

> This distribution is only temporarily available on a separate instance of the [Python Package index](https://test.pypi.org/project/sandy/).

## Running SANDY in a python shell
1. Open a terminal/Anaconda Prompt.
2. Switch to the correct python environment (if any).
3. Type `python` or `ipython` to open up a python/ipython session.
4. Then type
```python
import sandy
```

Now, if `sandy` is available you are ready to roll!

Still, many of the features of SANDY require running NJOY.
Here we assume that 


## Installation using Anaconda

### Downloading Anaconda
The recommended way to install the SANDY requirements (especially on Windows) is with the Anaconda distribution, which already includes several python packages and few other cool tools. The installer for the most recent Anaconda distributions can be found [here](https://www.anaconda.com/products/individual). For the Anaconda installation documentation, click [here](https://docs.anaconda.com/anaconda/install/).

### Getting started with Anaconda
The user guide to get yourself acquainted with Anaconda is available [here](https://docs.anaconda.com/anaconda/user-guide/). On a Windows machine we recommend that you use the Anaconda Prompt terminal.

# Installation from source
The source code of SANDY is available as a git repository. The recommended way to get SANDY is by cloning the source package using your local git distribution (*click [here](https://git-scm.com/downloads) to download it!*).

```sh
git clone https://github.com/luca-fiorito-11/sandy.git
cd sandy
```

If a `SSL Certificate problem` appears when pushing or pulling with git, one can tell git to not perform the validation of the certificate using the global option:

```git
git config --global http.sslVerify false
```

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

Now you can simply add the `sandy` folder to your `PYTHONPATH` and you are ready to go.

**in Linux**
```sh
export PYTHONPATH="/path/to/sandy:${PYTHONPATH}"
```
> Add the previous command to your `${HOME}/.bashrc` or `${HOME}/.bash_profile` file to make your changes permanent.

**in Windows**
```dos
set PYTHONPATH=%PYTHONPATH%;C:\path\to\sandy
```
> Alternatively, use the environment variable editor in **Control panel** &rarr; **View advanced system settings**.

Alternatively, you can add sandy to the python packages of your environent (although not necessary). Then, once you are in the folder `sandy`, run the installation command,
```sh
python setup.py install --user
```
> Don't forget to activate the virtual environment if you are using one.

SANDY is installed! Now go back to run it in a python shell.

## Testing SANDY
If you installed SANDY from source you can test if it works correctly by runnning a number (far away from being extensive enough) of unit tests. First, make sure you istalled the following packages

```sh
conda install -y --name sandy-devel -c conda-forge numpydoc nbval
```

Then, you can run the tests from a terminal/Anaconda Prompt from the `sandy` folder with
```sh
conda activate sandy-devel
mkdir tests
cd tests
pytest ../sandy
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
The installation of NJOY is not mandatory (it is if you want to produce random samples), but recommended.

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

If you run the following you should see `NJOY` in your list of environment variables (inside the python virtual environment),
```dos
conda env config vars list
```

To succesfully run NJOY Windows must be able to find some DLL files such as `cygwin1.dll`.
This file is part of cygwin, so most likely it's located in `C:\path\to\cygwin64\bin`.
Then, you have to add `C:\path\to\cygwin64\bin` (or the location where `cygwin1.dll` can be found) to your `PATH` typing the following on an Anaconda Prompt terminal 
```dos
set PATH=%PATH%;C:\path\to\cygwin64\bin
```
> Again, the environment variable editor in the control panel can also be used.


To verify that NJOY is correctly found by SANDY open a python terminal, import sandy and run `get_njoy()`.
```python
import sandy
sandy.get_njoy()
```

If you want to succesfully run NJOY2016 through SANDY, cygwin must be allowed to access different directories outside the cygwin home directory.
From a cygwin terminal this can be done specifing `/cygdrive/` before the directory absolute path.

> Example: you can access the root of your C: drive from cygwin by specifying the directory
```bash
cd /cygdrive/c
```

To consistently be able to write in your user account on the C: drive (`C:\Users\your_username`), for example, we recommend creating a symbolic link in a cygwin terminal, as
```bash
ln -sv /cygdrive/c/Users/your_username ~/your_username
```

To check if NJOY works you can:
 * open a Anaconda Prompt terminal;
 * move to a directory where cygwin has writing permission;
 * make sure `C:\path\to\cygwin64\bin` is in your list of environemnt variables (juste type `set` and look at `Path`);
 * run `C:\path\to\cygwin64\home\your_username\NJOY2016\bin\njoy.exe`.

NJOY should now be running in your terminal!