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
* pytest
* matplotlib
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
sandy --version
```

#### <a name="running-the-tests"></a>Running the tests

Once the installation is completed, type the following from any directory to automatically run SANDY's automated tests

```bash
sandy_tests
```

```sandy_tests``` accepts any ```pytest``` command line option.
Therefore, add the following options to:
 - ```-s``` to disable all captures
 - ```-v``` to display more information
 - ```--durations=N``` to get a list of the slowest 10 test durations
 
 ```bash
# Example
sandy_tests -s -v --durations=10
```

#### <a name="usage"></a>Usage

For an overview of SANDY's usage type

```bash
sandy --help
```

The command output is reported below

```
usage: sandy [-h] (--endf6-cov ENDF6_COV | --errorr-cov ERRORR_COV)
             [--samples SAMPLES] [--outdir OUTDIR] [-np PROCESSES] [--eig N]
             [-e ENERGY_POINT] [-v]
             file

Run SANDY

positional arguments:
  file                  ENDF-6 or PENDF format file.

optional arguments:
  -h, --help            show this help message and exit
  --endf6-cov ENDF6_COV
                        ENDF-6 file containing covariances.
  --errorr-cov ERRORR_COV
                        ERRORR file containing covariances.
  --samples SAMPLES     Number of samples.
  --outdir OUTDIR       Target directory where outputs are stored (default =
                        current working directory). If it does not exist it
                        will be created.
  -np PROCESSES, --processes PROCESSES
                        Number of worker processes (default=1).
  --eig N               Print the first N eigenvalues of the evaluated
                        covariance matrices (default = 0, do not print).
  -e ENERGY_POINT, --energy-point ENERGY_POINT
                        Additional energy points (in eV) to include in the
                        incoming-neutron energy grid (default = None). Provide
                        each energy point as an individual optional argument,
                        e.g. -e 100.0 -e 201.5
  -v, --version         SANDY's version.
```

#### <a name="contacts"></a>Contacts

* [**Luca Fiorito**](https://github.com/luca-fiorito-11) - lucafiorito.11@gmail.com

#### <a name="acknowledgemnts"></a>Acknowledgments

* SCKCEN
* ULB

#### <a name="publications"></a>Publications








