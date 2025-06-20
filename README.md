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
 - energy distrbutions of outgoing particles;
 - fission neutron multiplicities;
 - fission yields;
 - radioactive decay data.

### API for ENDF-6 files
The recent development on SANDY extended the original goal and focused on providing a simple interface for nuclear data files in ENDF-6 format.
Nuclear data such as cross sections, fission yields, radioactive decay constants and so on can be imported into tabulated dataframes (making extensive use of `pandas`) for further post-processing, analysis, plotting, ...

Examples are available [here](https://luca-fiorito-11.github.io/sandy_notebooks/).
 
 ***

## :wrench: Installation

SANDY can be installed both on Linux (recommended) or Windows.
The installation instructions are available [here](https://github.com/luca-fiorito-11/sandy/blob/develop/INSTALL.md).

<br>

## :hourglass: Development history and releases

The [latest](https://github.com/luca-fiorito-11/sandy/releases/latest) and older releases of SANDY are available [here](https://github.com/luca-fiorito-11/sandy/releases). 

For a detailed list of changes across versions, please refer to the [CHANGELOG](https://github.com/luca-fiorito-11/sandy/blob/develop/CHANGELOG.md) file.

<br>

## :notebook_with_decorative_cover: Documentation and how to use SANDY

The official SANDY documentation can be found [here](https://luca-fiorito-11.github.io/sandy-docs/index.html).

The primary use for SANDY is to produce perturbed nuclear data files that statistically represent the covariance information found in evaluated libraries.
This can be done using a command line interface.

### Example for cross sections, nubar and pfns sampling

```bash
# Produce ENDF-6 / PENDF perturbed files
python -m sandy.sampling  U235.jeff33  --processes 20  --samples 200

# Produce ACE perturbed files
python -m sandy.sampling  U235.jeff33  --processes 20  --samples 200  --acer  --temperatures 900 
```

For a more advanced use, look at these notebooks:
- [Direct sampling procedure: produce random ENDF-6 / ACE files](https://nbviewer.org/github/luca-fiorito-11/sandy_notebooks/blob/executed_notebooks_v1.1/notebook_random_files.ipynb)
- [2-step sampling procedure](https://github.com/luca-fiorito-11/sandy_notebooks/blob/executed_notebooks_v1.1/notebook_random_files_2steps.ipynb)

### Example for radioactive decay data sampling

```bash
python -m sandy.sampling  decay_data.jeff33  --processes 20  --samples 200
```

### Example for fission yield data sampling

```bash
# Only variance
python -m sandy.sampling nfy.jeff33 --processes 20 --samples 200

# With CEA covariance matrices for U235th and Pu239th 
python -m sandy.sampling nfy.jeff33 --processes 20 --samples 200  --fycov
```

For a more advanced use, look at these notebooks:
- [FY Sampling: create perturbed FY files (only variance, all energies)](https://nbviewer.org/github/luca-fiorito-11/sandy_notebooks/blob/executed_notebooks_v1.1/notebook_sampling_ify.ipynb)
- [FY sampling: create perturbed FY files with ad-hoc covariance matrix (all energies)](https://nbviewer.org/github/luca-fiorito-11/sandy_notebooks/blob/executed_notebooks_v1.1/notebook_sampling_ify_with_cov.ipynb)
- [FY sampling: create perturbed FY files using CEA conservative evaluation for U-235 thermal fission](https://nbviewer.org/github/luca-fiorito-11/sandy_notebooks/blob/executed_notebooks_v1.1/notebook_sampling_ify_u235th_cea_c1.ipynb)
- [FY sampling: create perturbed FY files using CEA conservative evaluation for Pu-239 thermal fission](https://nbviewer.org/github/luca-fiorito-11/sandy_notebooks/blob/executed_notebooks_v1.1/notebook_sampling_ify_pu239th_cea.ipynb)

<br>

## :video_game: Jupyter notebooks

[Here](https://luca-fiorito-11.github.io/sandy_notebooks/) you can find some cool [Jupyter notebooks](https://jupyter.org/) that kind of give an idea of what one can do with SANDY.

<br>

## :telephone_receiver: Contacts

[**Luca Fiorito**](https://github.com/luca-fiorito-11)

*  lucafiorito.11@gmail.com
*  luca.fiorito@sckcen.be

<br>

## :bookmark: Acknowledgments

SANDY was conceived and developed as a part of the PhD thesis on *Nuclear data uncertainty propagation and uncertainty quantification in nuclear codes* in the framework of a collaboration between [SCK CEN](https://www.sckcen.be) and [ULB](http://www.ulb.ac.be).

<br>

## :clipboard: Reference

Among the publications about SANDY, please use the following as references for citation.

### Burnup analysis
 - L. Fiorito, L. Engelen, F. Grimaldi, P. Romojaro, [*Nuclear data uncertainty propagation in the ARIANE GU3 burnup model using SANDY: Comparison between a FA and a pincell model*](https://doi.org/10.1016/j.anucene.2025.111423), Annals of Nuclear Energy, Volume 218, 2025, 111423, ISSN 0306-4549.

### Criticality
 - L. Fiorito, J. Dyrda and M. Fleming, [*JEFF-3.3 covariance application to ICSBEP using SANDY and NDAST*](https://doi.org/10.1051/epjconf/201921107003), EPJ Web of Conferences 211, 07003 (2019)

### Original publication
 - L. Fiorito, G. Žerovnik, A. Stankovskiy, G. Van den Eynde, P.E. Labeau, [*Nuclear data uncertainty propagation to integral responses using SANDY*](http://www.sciencedirect.com/science/article/pii/S0306454916305278), Annals of Nuclear Energy, Volume 101, 2017, Pages 359-366, ISSN 0306-4549.


<br>

## :earth_africa: Publications

This is a (incomplete) list of scientific studies citing SANDY.

- **Ebiwonjumi, Bamidele**. *Uncertainty Analyses of Tritium Production and Gamma Heating Rates in*. FUSION SCIENCE AND TECHNOLOGY, 2025. Article; Early Access, Vol. N/A. [DOI](https://doi.org/10.1080/15361055.2025.2498229)
- **Delipei, Gregory K.**. *Uncertainty Quantification Framework for High-Temperature Gas-Cooled*. NUCLEAR SCIENCE AND ENGINEERING, 2025. Article; Early Access, Vol. N/A. [DOI](https://doi.org/10.1080/00295639.2025.2486901)
- **Yaseen, Mahmoud**. *Sensitivity and uncertainty analysis in pebble-bed reactors: A study*. ANNALS OF NUCLEAR ENERGY, 2025. Article, Vol. 219. [DOI](https://doi.org/10.1016/j.anucene.2025.111428)
- **Fiorito, Luca**. *Nuclear data uncertainty propagation in the ARIANE GU3 burnup model*. ANNALS OF NUCLEAR ENERGY, 2025. Article, Vol. 218. [DOI](https://doi.org/10.1016/j.anucene.2025.111423)
- **Ryzhkov, Alexander A.**. *A review of the current nuclear data performance assessments in advanced*. ANNALS OF NUCLEAR ENERGY, 2025. Review, Vol. 212. [DOI](https://doi.org/10.1016/j.anucene.2024.110806)
- **Lovecky, M.**. *OPOS-1000: Advancing the efficiency of VVER-1000 spent nuclear fuel cask*. NUCLEAR ENGINEERING AND DESIGN, 2025. Article, Vol. 431. [DOI](https://doi.org/10.1016/j.nucengdes.2024.113723)
- **Dagan, R.**. *Investigation of nuclide inventory of cladding material irradiated in*. ANNALS OF NUCLEAR ENERGY, 2025. Article, Vol. 212. [DOI](https://doi.org/10.1016/j.anucene.2024.111061)
- **Lovecky, M.**. *Optimizing spent nuclear fuel cask loading for VVER-440 fuel*. NUCLEAR ENGINEERING AND TECHNOLOGY, 2024. Article, Vol. 56. [DOI](https://doi.org/10.1016/j.net.2024.07.014)
- **Osman, W.**. *AN INTEGRATED FRAMEWORK FOR UNCERTAINTY QUANTIFICATION IN HIGH*. Arxiv, 2024. preprint, Vol. N/A. [DOI](https://doi.org/arXiv:2411.03329)
- **Jo, YuGwon**. *Uncertainty quantification based on similarity analysis of reactor*. NUCLEAR ENGINEERING AND TECHNOLOGY, 2024. Article, Vol. 56. [DOI](https://doi.org/10.1016/j.net.2024.04.014)
- **Zu, Tiejun**. *Development and Verification of Sampling Code NECP-SOUL for Evaluated*. NUCLEAR SCIENCE AND ENGINEERING, 2025. Article, Vol. 199. [DOI](https://doi.org/10.1080/00295639.2024.2364471)
- **Fiorito, Luca**. *Uncertainty Quantification for the Doppler Reactivity Feedback*. NUCLEAR SCIENCE AND ENGINEERING, 2025. Review, Vol. 199. [DOI](https://doi.org/10.1080/00295639.2024.2353987)
- **Belfiore, Enrica**. *On the Assumptions Behind Statistical Sampling: A <SUP>235</SUP>U*. NUCLEAR SCIENCE AND ENGINEERING, 2025. Article, Vol. 199. [DOI](https://doi.org/10.1080/00295639.2024.2323217)
- **Lovecky, Martin**. *TOTAL MONTE CARLO UNCERTAINTY ANALYSIS OF VVER-440 SPENT NUCLEAR FUEL IN*. PROCEEDINGS OF 2024 31ST INTERNATIONAL CONFERENCE ON NUCLEAR, 2024. Proceedings Paper, Vol. N/A. [DOI](#)
- **Belanger, Hunter**. *Comparison and evaluation of resolved resonance region covariance*. JOURNAL OF NUCLEAR SCIENCE AND TECHNOLOGY, 2024. Article, Vol. 61. [DOI](https://doi.org/10.1080/00223131.2023.2279303)
- **Schnabel, Georg**. *How to explain ENDF-6 to computers: A formal ENDF format description*. Arxiv, 2023. preprint, Vol. N/A. [DOI](https://doi.org/arXiv:2312.08249)
- **Tada, Kenichi**. *Development of nuclear data processing code FRENDY version 2*. JOURNAL OF NUCLEAR SCIENCE AND TECHNOLOGY, 2024. Article, Vol. 61. [DOI](https://doi.org/10.1080/00223131.2023.2278600)
- **Lovecky, M.**. *Radiation damage analysis of the first generation VVER spent nuclear*. ANNALS OF NUCLEAR ENERGY, 2024. Article, Vol. 196. [DOI](https://doi.org/10.1016/j.anucene.2023.110214)
- **Martin-Hernandez, Guido**. *Device and method for low-uncertainty and high-efficiency neutron*. NUCLEAR INSTRUMENTS & METHODS IN PHYSICS RESEARCH SECTION A-ACCELERATORS, 2023. Article, Vol. 1057. [DOI](https://doi.org/10.1016/j.nima.2023.168712)
- **Lindley, B. A.**. *Ranking of nuclear data contributions to uncertainties in core physics*. ANNALS OF NUCLEAR ENERGY, 2023. Article, Vol. 194. [DOI](https://doi.org/10.1016/j.anucene.2023.110111)
- **Kleedtke, Noah**. *Utilization of ACE nuclear data file toolkit ACEtk to calculate relative*. ANNALS OF NUCLEAR ENERGY, 2023. Article, Vol. 193. [DOI](https://doi.org/10.1016/j.anucene.2023.110031)
- **Abrate, Nicolo**. *Nuclear Data Uncertainty Propagation for the Molten Salt Fast Reactor*. NUCLEAR SCIENCE AND ENGINEERING, 2023. Article, Vol. 197. [DOI](https://doi.org/10.1080/00295639.2023.2190861)
- **Grimaldi, Federico**. *Nuclear data uncertainty quantification on PWR spent nuclear fuel as a*. FRONTIERS IN ENERGY RESEARCH, 2023. Article, Vol. 11. [DOI](https://doi.org/10.3389/fenrg.2023.1146598)
- **Aimetta, Alex**. *A Nonintrusive Nuclear Data Uncertainty Propagation Study for the ARC*. NUCLEAR SCIENCE AND ENGINEERING, 2023. Article, Vol. 197. [DOI](https://doi.org/10.1080/00295639.2022.2153638)
- **Tada, Kenichi**. *Development of ACE file perturbation tool using FRENDY*. JOURNAL OF NUCLEAR SCIENCE AND TECHNOLOGY, 2023. Article, Vol. 60. [DOI](https://doi.org/10.1080/00223131.2022.2130463)
- **Ebiwonjumi, Bamidele**. *Propagation of radiation source uncertainties in spent fuel cask*. NUCLEAR ENGINEERING AND TECHNOLOGY, 2022. Article, Vol. 54. [DOI](https://doi.org/10.1016/j.net.2022.03.001)
- **Salino, Vivian**. *Incertitudes et ajustements de données nucléaires au moyen de méthodes*. , 2022. Dissertation/Thesis, Vol. N/A. [DOI](#)
- **Yamano, Naoki**. *Crucial importance of correlation between cross sections and angular*. JOURNAL OF NUCLEAR SCIENCE AND TECHNOLOGY, 2022. Article, Vol. 59. [DOI](https://doi.org/10.1080/00223131.2021.1997665)
- **Neudecker, Denise**. *Informing nuclear physics via machine learning methods with differential*. PHYSICAL REVIEW C, 2021. Article, Vol. 104. [DOI](https://doi.org/10.1103/PhysRevC.104.034611)
- **Patel, Vishal**. *An uncertainty quantification method relevant to material test reactors*. ANNALS OF NUCLEAR ENERGY, 2022. Article, Vol. 165. [DOI](https://doi.org/10.1016/j.anucene.2021.108629)
- **Abrate, Nicolo**. *Generalized perturbation techniques for uncertainty quantification in*. ANNALS OF NUCLEAR ENERGY, 2021. Article, Vol. 164. [DOI](https://doi.org/10.1016/j.anucene.2021.108623)
- **Gray, Ander**. *Uncertainty Propagation in SINBAD Fusion Benchmarks with Total Monte*. FUSION SCIENCE AND TECHNOLOGY, 2021. Article, Vol. 77. [DOI](https://doi.org/10.1080/15361055.2021.1895667)
- **Kos, Bor**. *ASUSD nuclear data sensitivity and uncertainty program package:*. NUCLEAR ENGINEERING AND TECHNOLOGY, 2021. Article, Vol. 53. [DOI](https://doi.org/10.1016/j.net.2021.01.034)
- **Skarbeli, Aris V.**. *Comparison of nuclear data uncertainties with other nuclear fuel cycle*. INTERNATIONAL JOURNAL OF ENERGY RESEARCH, 2021. Article, Vol. 45. [DOI](https://doi.org/10.1002/er.6992)
- **Fiorito, L.**. *On the use of criticality and depletion benchmarks for verification of*. ANNALS OF NUCLEAR ENERGY, 2021. Article, Vol. 161. [DOI](https://doi.org/10.1016/j.anucene.2021.108415)
- **Park, Jin Hun**. *Statistical Analysis of Tritium Breeding Ratio Deviations in the DEMO*. APPLIED SCIENCES-BASEL, 2021. Article, Vol. 11. [DOI](https://doi.org/10.3390/app11115234)
- **Neudecker, Denise**. *Which nuclear data can be validated with LLNL pulsed-sphere experiments?*. ANNALS OF NUCLEAR ENERGY, 2021. Article, Vol. 159. [DOI](https://doi.org/10.1016/j.anucene.2021.108345)
- **Skarbeli, Aris, V**. *Uncertainty and Optimization Analysis of Advanced Nuclear Fuel Cycles*. 30TH INTERNATIONAL CONFERENCE NUCLEAR ENERGY FOR NEW EUROPE (NENE 2021), 2021. Proceedings Paper, Vol. N/A. [DOI](#)
- **Fleming, M.**. *The High-Energy Intra-Nuclear Cascade Liege-based Residual (HEIR)*. ND 2019: INTERNATIONAL CONFERENCE ON NUCLEAR DATA FOR SCIENCE AND, 2020. Proceedings Paper, Vol. 239. [DOI](https://doi.org/10.1051/epjconf/202023920001)
- **Fleming, M.**. *New features and improvements in the NEA nuclear data tool suite*. ND 2019: INTERNATIONAL CONFERENCE ON NUCLEAR DATA FOR SCIENCE AND, 2020. Proceedings Paper, Vol. 239. [DOI](https://doi.org/10.1051/epjconf/202023919003)
- **Solis, Augusto Hernandez**. *Depletion uncertainty analysis to the MYRRHA fuel assembly model*. ND 2019: INTERNATIONAL CONFERENCE ON NUCLEAR DATA FOR SCIENCE AND, 2020. Proceedings Paper, Vol. 239. [DOI](https://doi.org/10.1051/epjconf/202023912001)
- **Sui, Zhuojie**. *Covariance-oriented sample transformation: A new sampling method for*. ANNALS OF NUCLEAR ENERGY, 2019. Article, Vol. 134. [DOI](https://doi.org/10.1016/j.anucene.2019.07.001)
- **Romojaro, P.**. *Sensitivity methods for effective delayed neutron fraction and neutron*. ANNALS OF NUCLEAR ENERGY, 2019. Article; Proceedings Paper, Vol. 126. [DOI](https://doi.org/10.1016/j.anucene.2018.11.042)
- **Goricanec, Tanja**. *Evaluation of the criticality and reaction rate benchmark experiments*. PROGRESS IN NUCLEAR ENERGY, 2019. Article, Vol. 111. [DOI](https://doi.org/10.1016/j.pnucene.2018.10.024)
- **Cabellos, Oscar**. *Examples of Monte Carlo techniques applied for nuclear data uncertainty*. 5TH INTERNATIONAL WORKSHOP ON NUCLEAR DATA EVALUATION FOR REACTOR, 2019. Proceedings Paper, Vol. 211. [DOI](https://doi.org/10.1051/epjconf/201921107008)
- **Fiorito, Luca**. *JEFF-3.3 covariance application to ICSBEP using SANDY and NDAST*. 5TH INTERNATIONAL WORKSHOP ON NUCLEAR DATA EVALUATION FOR REACTOR, 2019. Proceedings Paper, Vol. 211. [DOI](https://doi.org/10.1051/epjconf/201921107003)
- **Sjostrand, Henrik**. *Monte Carlo integral adjustment of nuclear data libraries - experimental*. 5TH INTERNATIONAL WORKSHOP ON NUCLEAR DATA EVALUATION FOR REACTOR, 2019. Proceedings Paper, Vol. 211. [DOI](https://doi.org/10.1051/epjconf/201921107007)
- **Kos, Bor**. *Validation and Use of Coupling SUSD3D with Denovo for Complex*. 28TH INTERNATIONAL CONFERENCE NUCLEAR ENERGY FOR NEW EUROPE (NENE 2019), 2019. Proceedings Paper, Vol. N/A. [DOI](#)
- **Stankovskiy, Alexey**. *High-energy nuclear data uncertainties propagated to MYRRHA safety*. ANNALS OF NUCLEAR ENERGY, 2018. Article, Vol. 120. [DOI](https://doi.org/10.1016/j.anucene.2018.05.041)
- **Trottier, Alexandre**. *Nuclear data sensitivity for reactor physics parameters in a lead-cooled*. ANNALS OF NUCLEAR ENERGY, 2018. Article, Vol. 120. [DOI](https://doi.org/10.1016/j.anucene.2018.05.047)
- **Castro, E.**. *Impact of the homogenization level, nodal or pin-by-pin, on the*. PROGRESS IN NUCLEAR ENERGY, 2018. Article, Vol. 104. [DOI](https://doi.org/10.1016/j.pnucene.2017.10.001)
- **Iwamoto, Hiroki**. *Monte Carlo uncertainty quantification of the effective delayed neutron*. JOURNAL OF NUCLEAR SCIENCE AND TECHNOLOGY, 2018. Article, Vol. 55. [DOI](https://doi.org/10.1080/00223131.2017.1416691)
- **Calic, Dusan**. *Nuclear Data Uncertainties of the BEAVRS Benchmark Core*. 27TH INTERNATIONAL CONFERENCE NUCLEAR ENERGY FOR NEW EUROPE (NENE 2018), 2018. Proceedings Paper, Vol. N/A. [DOI](#)
- **Kos, Bor**. *Coupling of the SUSD3D S/U Code With the Denovo Deterministic Transport*. 27TH INTERNATIONAL CONFERENCE NUCLEAR ENERGY FOR NEW EUROPE (NENE 2018), 2018. Proceedings Paper, Vol. N/A. [DOI](#)
- **Ambrozic, K.**. *Computational analysis of the dose rates at JSI TRIGA reactor*. APPLIED RADIATION AND ISOTOPES, 2017. Article, Vol. 130. [DOI](https://doi.org/10.1016/j.apradiso.2017.09.022)
- **Griseri, M.**. *Nuclear data uncertainty propagation on a sodium fast reactor*. NUCLEAR ENGINEERING AND DESIGN, 2017. Article, Vol. 324. [DOI](https://doi.org/10.1016/j.nucengdes.2017.08.018)
- **Stankovskiy, Alexey**. *Impact of intermediate and high energy nuclear data on the neutronic*. ND 2016: INTERNATIONAL CONFERENCE ON NUCLEAR DATA FOR SCIENCE AND, 2017. Proceedings Paper, Vol. 146. [DOI](https://doi.org/10.1051/epjconf/201714609001)
