[![Build Status](https://travis-ci.org/luca-fiorito-11/sandy.svg?branch=devel)](https://travis-ci.org/luca-fiorito-11/sandy)

# SANDY

SANDY is a python package that can read, write and perform a set of operations on nuclear data file in 
[ENDF-6 format](https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf).
The primay objective of the code, as it was originally conceived, is to produce *perturbed files* containing sampled parameters 
that represent the information stored in the evaluated nuclear data covariances.
Such files can be ultimately used to propagate uncertainties through any given compatible system using a brute force technique.

Currently, SANDY can draw samples for:
 - cross sections;
 - angular distrbutions of outgoing particles;
 - energy distrbutions of outgoing particles;
 - fission neutron multiplicities.
 
The sampling algorithm constructs multivariate normal distributions with a unit vector for mean and with relative 
covariances taken from the evaluated files.
Perturbation factors are sampled with the same multigroup structure of the covariance matrix, and are applied to the pointwise 
data to produce the perturbed files. 

## Content

- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Running the test](#running-the-tests)
- [Usage](#usage)
- [Contacts](#contacts)
- [Acknowledgments](#acknowledgments)
- [Publications](#publications)


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

SANDY is developed in `python3` and does not support `python2`.
In order to run SANDY, make sure that you have a version of python3 installed.

[Here](requirements.txt) you can find the python dependencies required to ensure the correct functioning of SANDY.

### Installation

To download SANDY, move to the folder where you want the source code and type

```git
git clone https://github.com/luca-fiorito-11/sandy.git
```

To install SANDY, run the following commands

```bash
cd sandy
python setup.py install
```

Make sure that `python` points to a correct `python3` executable for which you have administration rights.

To quickly check if SANDY was installed correctly, type the following from any directory

```bash
sandy --version
```

## Running the tests

Once the installation is completed, run ```pytest``` to automatically start SANDY's automated tests

```bash
pytest
```

More [pytest](https://docs.pytest.org/en/latest/) command line option can be added to customize the tests set.

## Usage

For an overview of SANDY's usage type

```bash
sandy --help
```


## <a name="contacts"></a>Contacts

* [**Luca Fiorito**](https://github.com/luca-fiorito-11) - lucafiorito.11@gmail.com

## Acknowledgments

SANDY was conceived and developed as a part of the PhD thesis on *Nuclear data uncertainty propagation and uncertainty quantification in nuclear codes* in the framework of a collaboration between [SCK-CEN](https://www.sckcen.be) and [ULB](http://www.ulb.ac.be).


#### <a name="publications"></a>Publications
L. Fiorito, G. Žerovnik, A. Stankovskiy, G. Van den Eynde, P.E. Labeau, [*Nuclear data uncertainty propagation to integral responses using SANDY*](http://www.sciencedirect.com/science/article/pii/S0306454916305278), Annals of Nuclear Energy, Volume 101, 2017, Pages 359-366, ISSN 0306-4549.






