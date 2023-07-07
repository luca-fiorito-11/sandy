<p align="center">
  <img src="./badges/python.svg" alt="Python version">
  <a href="https://travis-ci.org/luca-fiorito-11/sandy">
    <img src="https://travis-ci.org/luca-fiorito-11/sandy.svg?branch=master" alt="Build status">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
  <!---
  <a href="http://hits.dwyl.io/luca-fiorito-11/sandy">
    <img src="http://hits.dwyl.io/luca-fiorito-11/sandy.svg" alt="HitCount">
  </a>
-->
  <a href="https://coveralls.io/github/luca-fiorito-11/sandy">
    <img src="https://coveralls.io/repos/github/luca-fiorito-11/sandy/badge.svg" alt="Coverage Status">
  </a>
</p>


<h1 align="center" style="font-family:simplifica">SANDY</h1>


<h5 align="center">Sampling tool for nuclear data</h5>

<br>
 

SANDY is a python package that can read, write and perform a set of operations on nuclear data files in
[ENDF-6 format](https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf).

### Stochastic sampling of nuclear data 
The primary objective of the code, as it was originally conceived, is to produce *perturbed files* containing sampled parameters 
that represent the information stored in the evaluated nuclear data covariances.
Such files can be ultimately used to propagate uncertainties through any given compatible system using a brute force technique.

Currently, SANDY can draw samples for:
 - cross sections;
 - angular distrbutions of outgoing particles;
 - energy distrbutions of outgoing particles;
 - fission neutron multiplicities;
 - fission yields.

### API for ENDF-6 files
The recent development on SANDY extended the original goal and focused on providing a simple interface for nuclear data files in ENDF-6 format.
Nuclear data such as cross sections, fission yields, radioactive decay constants and so on can be imported into tabulated dataframes (making extensive use of `pandas`) for further post-processing, analysis, plotting, ...

Examples are available [here](https://luca-fiorito-11.github.io/sandy_notebooks/).
 
 ***

## :wrench: Installation

SANDY can be installed both on Linux (recommended) or Windows (using Anaconda).
The installation instructions are available [here](https://github.com/luca-fiorito-11/sandy/blob/develop/INSTALL.md).

<br>

## :hourglass: Development history and releases

The [latest](https://github.com/luca-fiorito-11/sandy/releases/latest) and older releases of SANDY are available [here](https://github.com/luca-fiorito-11/sandy/releases). 

<br>

## :notebook_with_decorative_cover: Documentation

The official SANDY documentation can be found [here](https://luca-fiorito-11.github.io/sandy-docs/index.html).

<br>

## :video_game: Jupyter notebooks

[Here](https://luca-fiorito-11.github.io/sandy_notebooks/) you can find some cool [Jupyter notebooks](https://jupyter.org/) that kind of give an idea of what one can do with SANDY. 

<br>

## :telephone_receiver: Contacts

* [**Luca Fiorito**](https://github.com/luca-fiorito-11) - lucafiorito.11@gmail.com

<br>

## :bookmark: Acknowledgments

SANDY was conceived and developed as a part of the PhD thesis on *Nuclear data uncertainty propagation and uncertainty quantification in nuclear codes* in the framework of a collaboration between [SCK CEN](https://www.sckcen.be) and [ULB](http://www.ulb.ac.be).

<br>

## :clipboard: Reference

Among the publications about SANDY, please use the following as references for citation.

 - L. Fiorito, J. Dyrda and M. Fleming, [*JEFF-3.3 covariance application to ICSBEP using SANDY and NDAST*](https://doi.org/10.1051/epjconf/201921107003), EPJ Web of Conferences 211, 07003 (2019)

 - L. Fiorito, G. Žerovnik, A. Stankovskiy, G. Van den Eynde, P.E. Labeau, [*Nuclear data uncertainty propagation to integral responses using SANDY*](http://www.sciencedirect.com/science/article/pii/S0306454916305278), Annals of Nuclear Energy, Volume 101, 2017, Pages 359-366, ISSN 0306-4549.

<br>

## :earth_africa: Publications

[Here](https://www.webofscience.com/wos/alldb/summary/4fbd3df7-2148-4510-b95c-91ac22b111b6-6bbbab39/date-descending/1) is a (incomplete) list of scientific studies citing SANDY.

> If some info are not correct or missing, please let us know!
