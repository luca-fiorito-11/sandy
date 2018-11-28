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
- [Reference](#reference)


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

> __Important__: SANDY works with files in ENDF-6 format. If you are not familiar with the format, have a look at the documentation [here](https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf).

### Prerequisites

SANDY is developed in `python3` and does not support `python2`.
In order to run SANDY, make sure that you have a version of python3 installed.

[Here](requirements.txt) you can find the python dependencies required to ensure the correct functioning of SANDY.

### Installation

To download SANDY, move to the folder where you want the source code and type

```bash
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

Once the installation is completed, run ```pytest``` to automatically start SANDY's tests

```bash
pytest
```

More [pytest](https://docs.pytest.org/en/latest/) command line option can be added to customize the tests set.

## Usage

For an overview of SANDY's usage type

```bash
sandy --help
```

## Examples


#### Data and covariances are in the same file

Produce 1000 perturbed copies of a ENDF-6 file `<tape>` that contains both evaluated data and covariances.
```
sandy  <tape>  --samples 1000
```

Below are reported the ENDF-6 data sections that will be perturbed and the respective covariance sections.

| Data type | Data section | Covariance section |
|----|:----:|:----:|
| fission multiplitcities | MF1 | MF31 |
| cross sections | MF3 | MF33 |
| angular ditributions | MF4 | MF34 |
| energy distributions | MF5 | MF35 |

> __Important__: cross sections will be perturbed __only__ if they are linearized and given in PENDF (pointwise-ENDF) format.
> To convert a ENDF-6 file into PENDF format, you can use nuclear data processing codes such as [NJOY](http://www.njoy21.io/NJOY2016/) or [PREPRO](https://www-nds.iaea.org/public/endf/prepro/).

#### Perturb only one or few data types

Add keyword option `--mf` to perturb only few data type.
For example, to produce 1000 perturbed copies of a file `<tape>` where only angular and energy distributions are perturbed, type
```
sandy  <tape>  --samples 1000  --mf 34 35
```

#### Data and covariances are in different files

Produce 1000 perturbed copies of a file `<tape>` that contains evaluated data using covariances from file `<covtape>`.
```
sandy  <tape>  --cov <covtape>  --samples 1000
```

> __Important__: this command is often used for perturbing cross sections, where the linearized data are in a PENDF file `<tape>` that might not contain covariances and the covariance data are in the original ENDF-6 file `<covtape>`.


#### Covariance data in ERRORR format

`ERRORR` is a [NJOY](http://www.njoy21.io/NJOY2016/) module that processes the covariance information present in a ENDF-6 file into a given multigroup structure. The resulting tabulated covariance is tabulated into an output file with a specific `ERRORR` format.
Not only does `ERRORR` process cross section covariances in MF33, but it can also handle the resonance-resonance covariances in ENDF-6 covariance section MF32.

To produce 1000 perturbed copies of a PENDF file `<tape>` including the MF32 covariances for resonance parameters, type
```
sandy  <tape>  --cov <covtape>  --samples 1000
```
where `<covtape>` is a `ERRORR` output file.



## <a name="contacts"></a>Contacts

* [**Luca Fiorito**](https://github.com/luca-fiorito-11) - lucafiorito.11@gmail.com

## Acknowledgments

SANDY was conceived and developed as a part of the PhD thesis on *Nuclear data uncertainty propagation and uncertainty quantification in nuclear codes* in the framework of a collaboration between [SCK-CEN](https://www.sckcen.be) and [ULB](http://www.ulb.ac.be).


## <a name="refrence"></a>Reference
Among the publications about SANDY, please use the following as a reference for citation.

L. Fiorito, G. Žerovnik, A. Stankovskiy, G. Van den Eynde, P.E. Labeau, [*Nuclear data uncertainty propagation to integral responses using SANDY*](http://www.sciencedirect.com/science/article/pii/S0306454916305278), Annals of Nuclear Energy, Volume 101, 2017, Pages 359-366, ISSN 0306-4549.






