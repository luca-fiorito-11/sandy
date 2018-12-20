![Alt text](./badges/python.svg)
[![Build Status](https://travis-ci.org/luca-fiorito-11/sandy.svg?branch=master)](https://travis-ci.org/luca-fiorito-11/sandy)
[![](https://img.shields.io/readthedocs/:packageName.svg)](https://luca-fiorito-11.github.io/sandy-docs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![HitCount](http://hits.dwyl.io/luca-fiorito-11/sandy.svg)](http://hits.dwyl.io/luca-fiorito-11/sandy)
[![Coverage Status](https://coveralls.io/repos/github/luca-fiorito-11/sandy/badge.svg)](https://coveralls.io/github/luca-fiorito-11/sandy)

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
 - fission neutron multiplicities;
 - fission yields.
 

## Installation

To install SANDY, run the following commands:

```
git clone https://github.com/luca-fiorito-11/sandy.git
cd sandy
python setup.py install
```

## Documentation

The official SANDY documentation can be found [here](https://luca-fiorito-11.github.io/sandy-docs/index.html).

## <a name="contacts"></a>Contacts

* [**Luca Fiorito**](https://github.com/luca-fiorito-11) - lucafiorito.11@gmail.com

## Acknowledgments

SANDY was conceived and developed as a part of the PhD thesis on *Nuclear data uncertainty propagation and uncertainty quantification in nuclear codes* in the framework of a collaboration between [SCK-CEN](https://www.sckcen.be) and [ULB](http://www.ulb.ac.be).


## <a name="refrence"></a>Reference
Among the publications about SANDY, please use the following as a reference for citation.

L. Fiorito, G. Žerovnik, A. Stankovskiy, G. Van den Eynde, P.E. Labeau, [*Nuclear data uncertainty propagation to integral responses using SANDY*](http://www.sciencedirect.com/science/article/pii/S0306454916305278), Annals of Nuclear Energy, Volume 101, 2017, Pages 359-366, ISSN 0306-4549.






