# Installation

## Download

The SANDY source code is available as a git repository.
The recommended way to get SANDY is by cloning the source package using your 
local git distribution 
(*click [here](https://git-scm.com/downloads) to download it!*).

```sh
git clone https://github.com/luca-fiorito-11/sandy.git
cd sandy
```

If a `SSL Certificate problem` appears when pushing or pulling with git, 
one can tell git to not perform the validation of the certificate using the 
global option:

```git
git config --global http.sslVerify false
```

## Installation using Anaconda

### Downloading Anaconda
The recommended way to install the SANDY requirements (especially on Windows) 
is with the Anaconda distribution, which already includes several python 
packages and few other cool tools.
The installer for the most recent Anaconda distributions can be found 
[here](https://www.anaconda.com/products/individual).
For the Anaconda installation documentation, click 
[here](https://docs.anaconda.com/anaconda/install/).

### Getting started with Anaconda
The user guide to get yourself acquainted with Anaconda is available 
[here](https://docs.anaconda.com/anaconda/user-guide/).
On a Windows machine we recommend that you use the Anaconda Prompt terminal.

### Installing SANDY

We advise to install the SANDY dependences in a python environment that was 
not previously altered.
You can do so from a terminal/Anaconda Prompt using the package/environment 
management system `conda` included in Anaconda, as
```sh
conda update --name base conda
conda create -y --name sandy-devel -c conda-forge python=3.7 jupyterlab jupyter_nbextensions_configurator
```
I like to include jupyterlab to my environment!

Then, dependences are installed from the requirement file,
```sh
conda install --force-reinstall -y --name sandy-devel -c conda-forge --file requirements_conda.txt
```
Make sure you are in the folder `sandy` that you cloned with git and where 
file `requirement.txt` is located.

From now on, every time you want to use SANDY in your python distribution 
you need to activate the correct environemt.
```sh
conda activate sandy-devel
...
conda deactivate
```
To manage your python environments read the 
[conda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).
Do not forget to add the `sandy` folder to your `PYTHONPATH`, or else you 
won't be able to use SANDY when you change folder.
You can do the following

**in Linux**
```sh
export PYTHONPATH="/path/to/sandy:${PYTHONPATH}"
```
Add the previous command to your `${HOME}/.bashrc` or `${HOME}/.bash_profile` 
file to make your changes permanent.

**in Windows**
```dos
set PYTHONPATH=%PYTHONPATH%;C:\path\to\sandy
```
or, if you want the changes to be permanent,
```dos
setx PYTHONPATH=%PYTHONPATH%;C:\path\to\sandy
```

## Running python
1. Open a terminal/Anaconda Prompt
2. Switch to the correct python environment
3. Type `python` or `ipython` to open up a python/ipython session,
4. Then type,
```python
import sandy
```

Now, if `sandy` is available you are ready to roll!

## Running SANDY in a Jupyter notebook

For combatibility issues we recommend installing a python kernel specific 
for the `sandy-devel` environment.
For that, you can run the following after making sure that `ipykernel` is 
installed in the virtual environment.

```sh
conda activate sandy-devel
python -m ipykernel install --user --name sandy-devel --display-name "Python3 (sandy-devel)"
```