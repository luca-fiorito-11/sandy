# Changelog

All notable changes to sandy are documented in this file.

## [1.1.0] - 2025-06-19
- SAdding Docker file, #378
- Modified parallelization for Windows compatibility, #365
- Parsing NJOY output for MF34, #359 
- Moved all modules at base level, #349, #350, 
- PFNS sampling (MF35), #335. #337, #338, #343, #356
- RDD sampling, #328, #329
- Modified way to filter by MT in covariance matrix without losing information, #327
- Adding URLs to retrieve online ND libraries, #321, #322
- Adding extra choices at command line interface level, #315, #367
- Adding energy-grids, #303
- FY sampling, including CEA covariance matrices for Pu239th and U235th , #293, #312, #314, #339, #341, #342, #344, #345, #352
- Read and write perturbation coefficients to and from Excel file for 2-step sampling (first `get_perturbations`, then `apply_perturbations`), #288 #290, #299
- Sampling from log-normal PDF as default, #289, #368
- Adding decay data files locally for faster parsing, #280
