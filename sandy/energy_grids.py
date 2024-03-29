import numpy as np

__all__ = [
    "CASMO12",
    "ECCO33",
    "ONEGROUP20",
    "SCALE238",
    "WIMS69",
    ]
__author__ = "Luca Fiorito"


ONEGROUP20 = np.array([
    1e-11,
    20,
    ]) * 1e6  # in eV


ECCO33 = np.array([
    1.00001e-05,
    1.00000e-01,
    5.40000e-01,
    4.00000e+00,
    8.31529e+00,
    1.37096e+01,
    2.26033e+01,
    4.01690e+01,
    6.79040e+01,
    9.16609e+01,
    1.48625e+02,
    3.04325e+02,
    4.53999e+02,
    7.48518e+02,
    1.23410e+03,
    2.03468e+03,
    3.35463e+03,
    5.53084e+03,
    9.11882e+03,
    1.50344e+04,
    2.47875e+04,
    4.08677e+04,
    6.73795e+04,
    1.11090e+05,
    1.83156e+05,
    3.01974e+05,
    4.97871e+05,
    8.20850e+05,
    1.35335e+06,
    2.23130e+06,
    3.67879e+06,
    6.06531e+06,
    1.00000e+07,
    1.96403e+07,
    ])   # in eV


WIMS69 = np.array([
    1e-11,
    5e-09,
    1e-08,
    1.5e-08,
    2e-08,
    2.5e-08,
    3e-08,
    3.5e-08,
    4.2e-08,
    5e-08,
    5.8e-08,
    6.7e-08,
    8e-08,
    1e-07,
    1.4e-07,
    1.8e-07,
    2.2e-07,
    2.5e-07,
    2.8e-07,
    3e-07,
    3.2e-07,
    3.5e-07,
    4e-07,
    5e-07,
    6.25e-07,
    7.8e-07,
    8.5e-07,
    9.1e-07,
    9.5e-07,
    9.72e-07,
    9.96e-07,
    1.02e-06,
    1.045e-06,
    1.071e-06,
    1.097e-06,
    1.123e-06,
    1.15e-06,
    1.3e-06,
    1.5e-06,
    2.1e-06,
    2.6e-06,
    3.3e-06,
    4e-06,
    9.877e-06,
    1.5968e-05,
    2.77e-05,
    4.8052e-05,
    7.5501e-05,
    0.00014873,
    0.00036726,
    0.0009069,
    0.0014251,
    0.0022395,
    0.0035191,
    0.00553,
    0.009118,
    0.01503,
    0.02478,
    0.04085,
    0.06734,
    0.111,
    0.183,
    0.3025,
    0.5,
    0.821,
    1.353,
    2.231,
    3.679,
    6.0655,
    10.0,
    ]) * 1e6  # in eV


CASMO12 = np.array([
    1.0000E-11,
    3.0000E-08,
    5.8000E-08,
    1.4000E-07,
    2.8000E-07,
    3.5000E-07,
    6.2500E-07,
    4.0000E-06,
    4.8052E-05,
    5.5300E-03,
    8.2100E-01,
    2.2310E+00,
    1.0000E+01,
    ]) * 1e6  # in eV


SCALE238 = np.array([
    1.00000E-11,
    1.00000E-10,
    5.00000E-10,
    7.50000E-10,
    1.00000E-09,
    1.20000E-09,
    1.50000E-09,
    2.00000E-09,
    2.50000E-09,
    3.00000E-09,
    4.00000E-09,
    5.00000E-09,
    7.50000E-09,
    1.00000E-08,
    2.53000E-08,
    3.00000E-08,
    4.00000E-08,
    5.00000E-08,
    6.00000E-08,
    7.00000E-08,
    8.00000E-08,
    9.00000E-08,
    1.00000E-07,
    1.25000E-07,
    1.50000E-07,
    1.75000E-07,
    2.00000E-07,
    2.25000E-07,
    2.50000E-07,
    2.75000E-07,
    3.00000E-07,
    3.25000E-07,
    3.50000E-07,
    3.75000E-07,
    4.00000E-07,
    4.50000E-07,
    5.00000E-07,
    5.50000E-07,
    6.00000E-07,
    6.25000E-07,
    6.50000E-07,
    7.00000E-07,
    7.50000E-07,
    8.00000E-07,
    8.50000E-07,
    9.00000E-07,
    9.25000E-07,
    9.50000E-07,
    9.75000E-07,
    1.00000E-06,
    1.01000E-06,
    1.02000E-06,
    1.03000E-06,
    1.04000E-06,
    1.05000E-06,
    1.06000E-06,
    1.07000E-06,
    1.08000E-06,
    1.09000E-06,
    1.10000E-06,
    1.11000E-06,
    1.12000E-06,
    1.13000E-06,
    1.14000E-06,
    1.15000E-06,
    1.17500E-06,
    1.20000E-06,
    1.22500E-06,
    1.25000E-06,
    1.30000E-06,
    1.35000E-06,
    1.40000E-06,
    1.45000E-06,
    1.50000E-06,
    1.59000E-06,
    1.68000E-06,
    1.77000E-06,
    1.86000E-06,
    1.94000E-06,
    2.00000E-06,
    2.12000E-06,
    2.21000E-06,
    2.30000E-06,
    2.38000E-06,
    2.47000E-06,
    2.57000E-06,
    2.67000E-06,
    2.77000E-06,
    2.87000E-06,
    2.97000E-06,
    3.00000E-06,
    3.05000E-06,
    3.15000E-06,
    3.50000E-06,
    3.73000E-06,
    4.00000E-06,
    4.75000E-06,
    5.00000E-06,
    5.40000E-06,
    6.00000E-06,
    6.25000E-06,
    6.50000E-06,
    6.75000E-06,
    7.00000E-06,
    7.15000E-06,
    8.10000E-06,
    9.10000E-06,
    1.00000E-05,
    1.15000E-05,
    1.19000E-05,
    1.29000E-05,
    1.37500E-05,
    1.44000E-05,
    1.51000E-05,
    1.60000E-05,
    1.70000E-05,
    1.85000E-05,
    1.90000E-05,
    2.00000E-05,
    2.10000E-05,
    2.25000E-05,
    2.50000E-05,
    2.75000E-05,
    3.00000E-05,
    3.12500E-05,
    3.17500E-05,
    3.32500E-05,
    3.37500E-05,
    3.46000E-05,
    3.55000E-05,
    3.70000E-05,
    3.80000E-05,
    3.91000E-05,
    3.96000E-05,
    4.10000E-05,
    4.24000E-05,
    4.40000E-05,
    4.52000E-05,
    4.70000E-05,
    4.83000E-05,
    4.92000E-05,
    5.06000E-05,
    5.20000E-05,
    5.34000E-05,
    5.90000E-05,
    6.10000E-05,
    6.50000E-05,
    6.75000E-05,
    7.20000E-05,
    7.60000E-05,
    8.00000E-05,
    8.20000E-05,
    9.00000E-05,
    1.00000E-04,
    1.08000E-04,
    1.15000E-04,
    1.19000E-04,
    1.22000E-04,
    1.86000E-04,
    1.92500E-04,
    2.07500E-04,
    2.10000E-04,
    2.40000E-04,
    2.85000E-04,
    3.05000E-04,
    5.50000E-04,
    6.70000E-04,
    6.83000E-04,
    9.50000E-04,
    1.15000E-03,
    1.50000E-03,
    1.55000E-03,
    1.80000E-03,
    2.20000E-03,
    2.29000E-03,
    2.58000E-03,
    3.00000E-03,
    3.74000E-03,
    3.90000E-03,
    6.00000E-03,
    8.03000E-03,
    9.50000E-03,
    1.30000E-02,
    1.70000E-02,
    2.50000E-02,
    3.00000E-02,
    4.50000E-02,
    5.00000E-02,
    5.20000E-02,
    6.00000E-02,
    7.30000E-02,
    7.50000E-02,
    8.20000E-02,
    8.50000E-02,
    1.00000E-01,
    1.28300E-01,
    1.50000E-01,
    2.00000E-01,
    2.70000E-01,
    3.30000E-01,
    4.00000E-01,
    4.20000E-01,
    4.40000E-01,
    4.70000E-01,
    4.99520E-01,
    5.50000E-01,
    5.73000E-01,
    6.00000E-01,
    6.70000E-01,
    6.79000E-01,
    7.50000E-01,
    8.20000E-01,
    8.61100E-01,
    8.75000E-01,
    9.00000E-01,
    9.20000E-01,
    1.01000E+00,
    1.10000E+00,
    1.20000E+00,
    1.25000E+00,
    1.31700E+00,
    1.35600E+00,
    1.40000E+00,
    1.50000E+00,
    1.85000E+00,
    2.35400E+00,
    2.47900E+00,
    3.00000E+00,
    4.30400E+00,
    4.80000E+00,
    6.43400E+00,
    8.18730E+00,
    1.00000E+01,
    1.28400E+01,
    1.38400E+01,
    1.45500E+01,
    1.56830E+01,
    1.73330E+01,
    2.00000E+01,
    ]) * 1e6  # in eV


ECCO33 = np.array([
    1.000010E-11,
    1.000000E-07,
    5.400000E-07,
    4.000000E-06,
    8.315287E-06,
    1.370959E-05,
    2.260329E-05,
    4.016900E-05,
    6.790405E-05,
    9.166088E-05,
    1.486254E-04,
    3.043248E-04,
    4.539993E-04,
    7.485183E-04,
    1.234098E-03,
    2.034684E-03,
    3.354626E-03,
    5.530844E-03,
    9.118820E-03,
    1.503439E-02,
    2.478752E-02,
    4.086771E-02,
    6.737947E-02,
    1.110900E-01,
    1.831564E-01,
    3.019738E-01,
    4.978707E-01,
    8.208500E-01,
    1.353353E+00,
    2.231302E+00,
    3.678794E+00,
    6.065307E+00,
    1.000000E+01,
    1.964033E+01,
    ]) * 1e6  # in eV
