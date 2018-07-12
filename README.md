[![Build Status](https://travis-ci.org/luca-fiorito-11/sandy.svg?branch=devel)](https://travis-ci.org/luca-fiorito-11/sandy)

# SANDY

### Sampling nuclear data and uncertainty

#### Installing SANDY

- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Running the test](#running-the-tests)
- [Usage](#usage)
- [Contacts](#contacts)
- [Acknowledgments](#acknowledgments)
- [Publications](#publications)


#### <a name="getting-started"></a>Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

##### <a name="prerequisites"></a>Prerequisites

SANDY is developed in `Python 3` and does not support `Python 2`.
In order to run SANDY, make sure that you have a version of `Python 3` installed.

The following dependencies are required to ensure the correct functioning of SANDY:

```
* pandas >= 0.20
* numpy
* scipy
* pytest >= 3.3
```

##### <a name="installation"></a>Installation

To download SANDY, move to the folder where you want the source code and type

```git
git clone https://github.com/luca-fiorito-11/sandy.git
```

To install SANDY, run the following commands

```bash
cd sandy
python setup.py install
```

Make sure that `python` points to a correct `Python 3` executable for which you have administration rights.

To quickly check if SANDY was installed correctly, type

```bash
python -m sandy.sampling --version
```

#### <a name="running-the-tests"></a>Running the tests

Once the installation is completed, run ```pytest``` to automatically start SANDY's automated tests

```bash
pytest
```

#### <a name="usage"></a>Usage

For an overview of SANDY's usage type

```bash
python -m sandy.sampling --help
```

The command output is reported below

```
usage: python -m sandy.sampling [-h]
                                (--endf6-cov ENDF6_COV | --errorr-cov ERRORR_COV)
                                --samples SAMPLES [--outdir DIR]
                                [-np PROCESSES] [--eig N]
                                [--mat {1,..,9999} [{1,..,9999} ...]]
                                [--mf {1,..,40} [{1,..,40} ...]]
                                [--mt {1,..,999} [{1,..,999} ...]]
                                [--outname OUTNAME] [--verbose] [-e E [E ...]]
                                [-v]
                                file

Run sampling

positional arguments:
  file                  ENDF-6 or PENDF format file

optional arguments:
  -h, --help            show this help message and exit
  --endf6-cov ENDF6_COV
                        ENDF-6 file containing covariances
  --errorr-cov ERRORR_COV
                        ERRORR file containing covariances
  --samples SAMPLES     number of samples
  --outdir DIR          target directory where outputs are stored
                        (default = current working directory)
                        if it does not exist it will be created
  -np PROCESSES, --processes PROCESSES
                        number of worker processes (default = 1)
  --eig N               print the first N eigenvalues of the evaluated covariance matrices
                        (default = do not print)
  --mat {1,..,9999} [{1,..,9999} ...]
                        draw samples only from the selected MAT sections (default = keep all)
  --mf {1,..,40} [{1,..,40} ...]
                        draw samples only from the selected MF sections (default = keep all)
  --mt {1,..,999} [{1,..,999} ...]
                        draw samples only from the selected MT sections (default = keep all)
  --outname OUTNAME     basename for the output files (default is the the basename of <file>.)
  --verbose             turn on verbosity (default = quiet)
  -e E [E ...], --energy-points E [E ...]
                        additional energy points (in eV) to include in the incoming-neutron energy grid
                        (default = None)
  -v, --version         SANDY's version.
```

#### <a name="contacts"></a>Contacts

* [**Luca Fiorito**](https://github.com/luca-fiorito-11) - lucafiorito.11@gmail.com

#### <a name="acknowledgemnts"></a>Acknowledgments

* SCKCEN
* ULB

#### <a name="publications"></a>Publications








