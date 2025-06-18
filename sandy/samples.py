import numpy as np
import pandas as pd
import os
from math import isnan
import logging
from numpy.linalg import norm as matrixnorm
from scipy.stats import kstest, norm

from .xs import Xs, redundant_xs

__author__ = "Luca Fiorito"
__all__ = [
        "Samples",
        "read_fy_samples",
        ]


def read_fy_samples(file='PERT_MF8_MT454.xlsx'):
    """
    Read relative perturbations for fission yields from excel file produced by
    :obj:`~sandy.endf6.Endf6.get_perturbations_fy`.

    Parameters
    ----------
    file : `str`, optional
        The name of the file containing the perturbations for MT454.
        The default is `'PERT_MF8_MT454.xlsx'`.
        The default of the tabulated excel file is:
            
            - 1st column: energy in eV
            - 2nd column: ZAP
            - 3rd - nth columns: sample ID
        
        The name of each sheet is a ZAM number.

    Returns
    -------
    smp : `pd.DataFrame`
        Dataframe with perturbation coefficients given per:
            
            - ZAM: fissioning nuclide
            - E: neutron energy
            - ZAP: fission product
            - SMP: sample ID

    Notes
    -----
    .. note:: This does not use the object :obj:`~sandy.samples.Samples`.


    Examples
    --------
    
    Default use case.
    Produce an excel file of samples (verbosity needed to produce the excel file).

    >>> import sandy
    >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", [922350, 922380])
    >>> smps = tape.get_perturbations(2, verbose=True)

    Read it.
    
    >>> smps2 = sandy.read_fy_samples()

    Test that it was read correctly.

    >>> smps = smps.astype({'ZAM': 'int32', 'E': 'float64', 'ZAP': 'int32', 'SMP': 'int32', 'VALS': 'float64'})
    >>> smps2 = smps2.astype({'ZAM': 'int32', 'E': 'float64', 'ZAP': 'int32', 'SMP': 'int32', 'VALS': 'float64'})
    >>> assert smps2[["ZAM", "E", "ZAP", "SMP"]].equals(smps[["ZAM", "E", "ZAP", "SMP"]])
    >>> np.testing.assert_array_almost_equal(smps2.VALS, smps.VALS)
    """
    all_sheets = pd.read_excel(file, sheet_name=None)
    
    smp = []
    for k, v in all_sheets.items():
        s = v.ffill().assign(ZAM=int(k)).set_index(["ZAM", "E", "ZAP"]).rename_axis("SMP", axis=1).stack().rename("VALS").reset_index()
        smp.append(s)
    
    # same sorting structure as when it was produced in get_perturbations_fy
    smp = pd.concat(smp, ignore_index=True).sort_values(by=["ZAM", "E", "ZAP", "SMP"])

    return smp



def summarize_sample(s, cov):
    """
    Compute statistical summaries of sample covariance deviations from a reference covariance matrix.

    This function quantifies deviations between the sample covariance matrix obtained 
    from `s` and the reference covariance matrix `cov` using matrix norms. It evaluates 
    the total deviation, diagonal deviation, and off-diagonal deviation. Additionally, 
    it assesses normality and log-normality of the sample distribution.

    Parameters
    ----------
    s : :obj:`~sandy.samples.Samples`
        An object containing sample data from which the sample covariance matrix is computed.
    cov : :obj:`~sandy.samples.CategoryCov`
        A DataFrame containing the reference covariance matrix for comparison.

    Returns
    -------
    dict
        A dictionary summarizing the deviations of the sample covariance matrix from 
        the reference covariance matrix across different metrics.

        - `"Sample Size"` : int  
          The number of samples used.
        - `"Frobenius Norm"` : float  
          Relative Frobenius norm of the deviation matrix.
        - `"Diag Frobenius Norm"` : float  
          Relative Frobenius norm of the deviation in the diagonal elements.
        - `"Off-Diag Frobenius Norm"` : float  
          Relative Frobenius norm of the deviation in the off-diagonal elements.
        - `"Frobenius Norm < 5%"` : bool  
          Whether the total deviation is below 5%.
        - `"Diag Frobenius Norm < 5%"` : bool  
          Whether the diagonal deviation is below 5%.
        - `"Off-Diag Frobenius Norm < 5%"` : bool  
          Whether the off-diagonal deviation is below 5%.
        - `"% Accept Normality Test"` : float  
          Percentage of rows passing the Kolmogorov-Smirnov normality test.
        - `"% Accept LogNormality Test"` : float  
          Percentage of rows passing the log-normality test.

    Notes
    -----
    - The function normalizes deviations using the Frobenius norm.
    - The diagonal and off-diagonal components are analyzed separately.
    - Normality and log-normality tests are performed on the sample data.

    Warnings
    --------
    If the sample size is less than 2, a warning is logged, and an empty dictionary is returned.

    """

    # Extract covariance matrices
    C = cov.data.values
    SC = s.get_cov().values
    N = s.data.shape[1]

    if N < 2:
        logging.warning(f"cannot produce sample summary because sample size is {N}")
        return {}

    # Compute reference norms
    ref_total = matrixnorm(C)
    ref_diag = matrixnorm(np.diag(np.diag(C)))
    ref_offdiag = matrixnorm(C - np.diag(np.diag(C)))
    
    diff_total = matrixnorm(SC - C)
    diff_diag = matrixnorm(np.diag(np.diag(SC - C)))
    diff_offdiag = matrixnorm((SC - C) - np.diag(np.diag(SC - C)))

    # test normality and lognomality
    test_norm = s.test_normality().values
    fraction_normal = test_norm.sum() / test_norm.size
    test_lognorm = s.test_normality(lognormal=True).values
    fraction_lognormal = test_lognorm.sum() / test_lognorm.size

    # Compute overall metrics
    summary = {
        "Sample Size": N,
        "Frobenius Norm Diff": diff_total / ref_total,
        "Diag Frobenius Norm Diff": diff_diag / ref_diag,
        "Off-Diag Frobenius Norm Diff": diff_offdiag / ref_offdiag,
        "Accept Normality Test": fraction_normal,
        "Accept LogNormality Test": fraction_lognormal,
        }
    return summary



class Samples():
    """
    A container for managing and analyzing Monte Carlo samples of nuclear data.

    This class provides a structured way to store, manipulate, and analyze samples 
    from nuclear data perturbations. It includes methods for computing statistical 
    properties such as mean, standard deviation, covariance, and correlation matrices.
    Additionally, it offers tools for data transformation, exporting, and normality 
    testing.

    Attributes
    ----------
    data : pd.DataFrame
        A DataFrame containing the sample values, where:
        - Rows represent individual nuclear reactions or energy bins.
        - Columns represent different sample realizations.

    Methods
    -------
    apply_function(func)
        Apply a function to each sample.
    from_excel(file, beg=None, end=None)
        Load sample data from an Excel file.
    get_corr()
        Compute and return the correlation matrix of the samples.
    get_cov()
        Compute and return the covariance matrix of the samples.
    get_eleft()
        Replace energy intervals with their left bounds.
    get_eright()
        Replace energy intervals with their right bounds.
    get_mean()
        Compute and return the mean values of the samples.
    get_std()
        Compute and return the standard deviation of the samples.
    get_rstd()
        Compute and return the relative standard deviation of the samples.
    iterate_xs_samples()
        Generator that iterates over each sample as an `Xs` object.
    test_normality(lognormal=False, alpha=0.05)
        Perform a Kolmogorov-Smirnov test for normality.
    to_excel(file)
        Export sample data to an Excel file.
    truncate_normal()
        Truncate negative and extreme values in normally distributed samples.
    """

    _columnsname = "SMP"

    def __init__(self, df, *args, **kwargs):
        self.data = pd.DataFrame(df, *args, dtype=float, **kwargs)

    def __repr__(self):
        return self.data.__repr__()

    @property
    def data(self):
        """
        Dataframe of samples.

        Attributes
        ----------
        index : `pd.Index` or `pandas.MultiIndex`
            indices
        columns : `pd.Index`
            samples numbering
        values : `np.array`
            sample values as `float`

        Returns
        -------
        `pd.DataFrame`
            tabulated samples
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data.rename_axis(self.__class__._columnsname, axis=1)
        self._data.columns = self._data.columns.astype(int)

    def get_mean(self):
        """
        Compute the mean of each row in the dataset.
    
        This method calculates the mean along axis 1 (row-wise) and returns a pandas Series.
    
        Returns
        -------
        pd.Series
            A Series containing the mean value for each row.
            The index matches `self.data.index`, and the series is named "MEAN".
    
        Examples
        --------
        >>> import sandy
        >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["row1", "row2"])
        >>> obj = sandy.Samples(df)
        >>> obj.get_mean()
        row1    2.00000e+00
        row2    5.00000e+00
        Name: MEAN, dtype: float64

        """
        return self.data.mean(axis=1).rename("MEAN")

    def get_corr(self):
        """
        Compute the correlation matrix of the dataset.
    
        This method calculates the Pearson correlation coefficients between columns 
        in `self.data` and returns a correlation matrix.
    
        Returns
        -------
        pd.DataFrame
            A square matrix where each entry (i, j) represents the Pearson correlation 
            coefficient between column i and column j.
    
        Examples
        --------
        >>> import sandy
        >>> df = pd.DataFrame({0: [1, 2, 1], 1: [2, 4, 3], 2: [3, 6, 2]}, index=["A", "B", "C"])
        >>> obj = sandy.Samples(df)
        >>> obj.get_corr()
                 A         B         C
        A  1.00000e+00  1.00000e+00  5.00000e-01
        B  1.00000e+00  1.00000e+00  5.00000e-01
        C  5.00000e-01  5.00000e-01  1.00000e+00
        """
        return self.data.T.corr()

    def get_cov(self):
        """
        Compute the covariance matrix of the dataset.
    
        This method calculates the covariance between columns in `self.data` and 
        returns a covariance matrix.
    
        Returns
        -------
        pd.DataFrame
            A square matrix where each entry (i, j) represents the covariance between 
            column i and column j.
    
        Examples
        --------
        >>> import sandy
        >>> df = pd.DataFrame({0: [1, 2, 1], 1: [2, 4, 3], 2: [3, 6, 2]}, index=["A", "B", "C"])
        >>> obj = sandy.Samples(df)
        >>> obj.get_cov()
                 A         B         C
        A  1.00000e+00  2.00000e+00  5.00000e-01
        B  2.00000e+00  4.00000e+00  1.00000e+00
        C  5.00000e-01  1.00000e+00  1.00000e+00
    
        """
        return self.data.T.cov()

    def get_eright(self):
        """
        Create a new DataFrame with an updated index, where the 'E' level in the MultiIndex 
        (representing energy intervals) is replaced by the right bound of each interval, 
        and the original 'E' level is dropped. The new level is named 'ERIGHT'.
    
        The method performs the following steps:
        1. Extracts the right bounds of the energy intervals in the 'E' level of the MultiIndex.
        2. Constructs a new MultiIndex that includes the right bound values as a new level, 
           and removes the original 'E' level.
        3. Creates a copy of the current DataFrame, replacing the index with the newly 
           constructed MultiIndex.
        4. Returns the modified DataFrame with the updated index.
    
        Returns
        -------
        pd.DataFrame
            A new DataFrame with the same data as the original, but with an updated MultiIndex 
            where the 'E' level has been replaced by the 'ERIGHT' level.
    
        Notes
        -----
        - The method assumes that the original MultiIndex has a level named 'E', which contains 
          intervals.
        - The method will drop the 'E' level from the MultiIndex and replace it with the right bounds 
          of the intervals, renaming the new level to 'ERIGHT'.

        Examples
        --------
        >>> import sandy
        >>> import numpy as np
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10010)
        >>> smps1 = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=2, ek=[1, 2, 3])))[33]
        >>> np.testing.assert_array_equal(smps1.get_eright().index.get_level_values("ERIGHT"), [2, 3])
        >>> np.testing.assert_array_equal(smps1.get_eright().values, smps1.data.values)

        """
        multi_index = self.data.index

        # right bound of the energy intervals, which will become the new level
        new_level = [x.right for x in multi_index.get_level_values("E")]

        # Add the new level and drop the EnergyInterval index
        new_multi_index = pd.MultiIndex.from_tuples(
            [(*old, new) for old, new in zip(multi_index, new_level)],
            names=[*multi_index.names, "ERIGHT"]
        ).droplevel("E")

        # Create a copy of the samples with the new index
        # deep=True prevents overwriting self
        copy = self.data.copy(deep=True)
        copy.index = new_multi_index

        # Return a dataframe, not the Samples object
        return copy

    def get_eleft(self):
        """
        Create a new DataFrame with an updated index, where the 'E' level in the MultiIndex 
        (representing energy intervals) is replaced by the left bound of each interval, 
        and the original 'E' level is dropped. The new level is named 'ELEFT'.
    
        The method performs the following steps:
        1. Extracts the left bounds of the energy intervals in the 'E' level of the MultiIndex.
        2. Constructs a new MultiIndex that includes the left bound values as a new level, 
           and removes the original 'E' level.
        3. Creates a copy of the current DataFrame, replacing the index with the newly 
           constructed MultiIndex.
        4. Returns the modified DataFrame with the updated index.
    
        Returns
        -------
        pd.DataFrame
            A new DataFrame with the same data as the original, but with an updated MultiIndex 
            where the 'E' level has been replaced by the 'ELEFT' level.
    
        Notes
        -----
        - The method assumes that the original MultiIndex has a level named 'E', which contains 
          intervals.
        - The method will drop the 'E' level from the MultiIndex and replace it with the left bounds 
          of the intervals, renaming the new level to 'ELEFT'.

        Examples
        --------
        >>> import sandy
        >>> import numpy as np
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10010)
        >>> smps1 = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=2, ek=[1, 2, 3])))[33]
        >>> np.testing.assert_array_equal(smps1.get_eleft().index.get_level_values("ELEFT"), [1, 2])
        >>> np.testing.assert_array_equal(smps1.get_eleft().values, smps1.data.values)

        """

        multi_index = self.data.index

        # right bound of the energy intervals, which will become the new level
        new_level = [x.left for x in multi_index.get_level_values("E")]

        # Add the new level and drop the EnergyInterval index
        new_multi_index = pd.MultiIndex.from_tuples(
            [(*old, new) for old, new in zip(multi_index, new_level)],
            names=[*multi_index.names, "ELEFT"]
        ).droplevel("E")

        # Create a copy of the samples with the new index
        # deep=True prevents overwriting self
        copy = self.data.copy(deep=True)
        copy.index = new_multi_index

        # Return a dataframe, not the Samples object
        return copy

    def get_std(self):
        """
        Compute the standard deviation for each row in the dataset.
    
        Returns
        -------
        pd.Series
            A Series containing the row-wise standard deviation of `self.data`, 
            named "STD". The index matches `self.data.index`.
    
        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> data = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=["A", "B", "C"])
        >>> obj = Samples(data)
        >>> obj.get_std()
        A   1.00000e+00
        B   1.00000e+00
        C   1.00000e+00
        Name: STD, dtype: float64

        """
        return self.data.std(axis=1).rename("STD")
    
    
    def get_rstd(self):
        """
        Compute the relative standard deviation (coefficient of variation) for each row.
    
        The relative standard deviation (RSTD) is calculated as the standard deviation 
        divided by the mean for each row. It quantifies variability relative to the mean.
    
        Returns
        -------
        pd.Series
            A Series containing the row-wise relative standard deviation, named "RSTD".
            The index matches `self.data.index`.
    
        Notes
        -----
        - If the mean of a row is zero, the result may contain NaN or infinite values.
    
        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> import sandy
        >>> data = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=["A", "B", "C"])
        >>> obj = sandy.Samples(data)
        >>> obj.get_rstd()
        A   5.00000e-01
        B   2.00000e-01
        C   1.25000e-01
        Name: RSTD, dtype: float64

        """
        return (self.get_std() / self.get_mean()).rename("RSTD")


    def apply_function(self, foo):
        """
        Apply a given function to all items in the sample.

        Parameters
        ----------
        foo: `func`
            Function to apply to each item.

        Returns
        -------
        :obj: `~sandy.samples.Samples`
            Samples with applied function.

        Examples
        --------
        
        Test that this method can be used to shift samples.

        >>> import sandy
        >>> smp = sandy.Samples([[1, 1], [1, 1]])
        >>> shifted_smp = smp.apply_function(lambda x: x + 1).data
        >>> np.testing.assert_array_equal(shifted_smp, np.ones((2, 2)) * 2)

        """
        samples = foo(self.data)
        return self.__class__(samples)

    def truncate_normal(self, mode=1):
        """
        This method assumes the samples are in relative units, hence centered in
        a vector of 1.

        Parameters
        ----------
        mode : `int`, optional
            truncation mode. The default is 1.
            `mode=1` sets values `<0` to `0`, and values `>2` to `2`.
            `mode=2` sets values `<0` and `>2` to `1`

        Raises
        ------
        ValueError
            Only mode=1 or mode=2 are allowed.

        Returns
        -------
        :obj: `~sandy.samples.Samples`
            Truncated samples.

        Examples
        --------
        
        Test mode 1.

        >>> import sandy
        >>> smp = sandy.Samples([[3, 1], [-2, 0.5]])
        >>> smp_1 = smp.truncate_normal(mode=1).data
        >>> np.testing.assert_array_equal(smp_1, [[2, 1], [0, 0.5]])

        Test mode 2.

        >>> smp_2 = smp.truncate_normal(mode=2).data
        >>> np.testing.assert_array_equal(smp_2, [[1, 1], [1, 0.5]])

        """
        samples = self.data

        if mode == 1:
            lower_bound = samples > 0
            upper_bound = samples < 2
            samples = samples.where(lower_bound, 0)
            samples = samples.where(upper_bound, 2)

        elif mode == 2:
            lower_bound = samples > 0
            upper_bound = samples < 2
            samples = samples.where(lower_bound, 1)
            samples = samples.where(upper_bound, 1)

        else:
            raise ValueError("only mode=1 or mode=2 are allowed")
        
        return self.__class__(samples)
    
    def test_normality(self, lognormal=False, alpha=0.05):
        """
        Perform a Kolmogorov-Smirnov (KS) test for normality on each row of the dataset.
    
        This method tests whether each row of the dataset follows a normal distribution 
        by standardizing the data and comparing it to a standard normal distribution. 
        Optionally, the data can be log-transformed before testing.
    
        Parameters
        ----------
        lognormal : bool, optional
            If True, applies a natural logarithm transformation to the data before testing.
            This is useful for detecting lognormal distributions. Default is False.
            
        alpha : float, optional
            Significance level for the KS test (default is 0.05).
            - If p-value > alpha: Fail to reject H0 (data appears normal).
            - If p-value ≤ alpha: Reject H0 (data significantly deviates from normality).
    
        Returns
        -------
        pd.Series
            A Series indicating whether each row passes the normality test (True = normal).
            The index matches the non-constant rows of `self.data.index`, and the series 
            is named "KS TEST".
    
        Notes
        -----
        - If `lognormal=True`, the natural logarithm of the data is used before testing.
        - Each row is standardized to have zero mean and unit variance before testing.
        - The KS test compares each row’s empirical distribution to a standard normal distribution.
        - Rows with constant values (zero standard deviation) are excluded from the test.
        - The result is a boolean series where `True` indicates normality.
    
        Examples
        --------

        >>> import numpy as np
        >>> import pandas as pd
        >>> from sandy.samples import Samples
        >>> rng = np.random.default_rng(42)
    
        Generate a normal dataset (should pass the normality test).
            
        >>> normal_data = rng.normal(loc=0, scale=1, size=(5, 1000))
        >>> samples = Samples(pd.DataFrame(normal_data))
        >>> normality_result = samples.test_normality()
        >>> assert normality_result.all()  # All rows should be normally distributed
    
        Generate a lognormal dataset (should fail normality test but pass lognormal test).

        >>> lognormal_data = np.exp(normal_data)
        >>> samples = Samples(pd.DataFrame(lognormal_data))
        >>> normality_result = samples.test_normality()
        >>> assert not normality_result.any()  # Lognormal data should fail the normality test
        >>> lognormality_result = samples.test_normality(lognormal=True)
        >>> assert lognormality_result.all()  # Lognormal data should pass when tested for lognormality

        """
        # Convert data to a NumPy array
        s = self.data.copy().values

        if lognormal:
            s = np.log(s)
        
        # Compute row-wise mean and standard deviation
        row_means = np.mean(s, axis=1, keepdims=True)
        row_stds = np.std(s, axis=1, keepdims=True)
        
        valid_rows = row_stds.squeeze() > 0  # Identify non-constant rows
        
        # Only standardize valid rows
        standardized_data = (s[valid_rows] - row_means[valid_rows]) / row_stds[valid_rows]
   
        # Perform KS test for each row and store p-values
        p_values = np.array([kstest(row, norm.cdf).pvalue for row in standardized_data])

        # Determine normality based on alpha threshold
        normality_results = p_values > alpha
    
        # Return results as a pandas Series
        series = pd.Series(normality_results, index=self.data.index[valid_rows], name="KS TEST")
        return series

    def to_excel(self, file):
        """
        Save the sample dataset to an Excel file.
    
        This method exports the sample data to an Excel file, adding separate columns 
        for the left and right bounds of the 'E' interval before saving. The data is 
        written to the 'SMP' sheet of the specified Excel file.
    
        Parameters
        ----------
        file : str
            The path to the Excel file. If the file exists, data is appended; otherwise, 
            a new file is created.
    
        Returns
        -------
        None
            The function does not return a value but writes data to an Excel file.
    
        Notes
        -----
        - The 'E' column, assumed to contain interval-like objects, is split into 
          'ELEFT' (left bound) and 'ERIGHT' (right bound) before being removed.
        - The data is saved in a sheet named 'SMP' using the 'openpyxl' engine.
        - If the file already exists, the function appends the data instead of overwriting it.

        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> from sandy.samples import Samples
        >>> from pandas import Interval
    
        Create sample data with an interval index.

        >>> data = pd.DataFrame(
        ...     {
        ...         "MT": [102, 102, 18, 18],
        ...         "E": [Interval(1e-5, 1.0), Interval(1.0, 20.0), Interval(1e-5, 1.0), Interval(1.0, 20.0)],
        ...         0: [1.01, 0.99, 1.02, 1.03],
        ...         1: [1.02, 1.00, 1.01, 1.05],
        ...     }
        ... ).set_index(["MT", "E"])
        
        >>> samples = Samples(data)
        >>> samples.to_excel("test_output.xlsx")
    
        Read the data back and check correctness.

        >>> df_read = pd.read_excel("test_output.xlsx", sheet_name="SMP")
        >>> assert "ELEFT" in df_read.columns and "ERIGHT" in df_read.columns
        >>> assert df_read.shape == (4, 5)  # Expect 4 rows and 5 columns including ELEFT and ERIGHT

        """
        # Reset index and prepare data
        # deep=True makes also a copy of the metadata, or else the index of self will be change
        df = self.data.copy(deep=True)
    
        if "E" in df.index.names:
            index_df = df.index.to_frame()

            # Replace "E" with "ELEFT" and "ERIGHT"
            index_df["ELEFT"] = index_df["E"].apply(lambda x: x.left)
            index_df["ERIGHT"] = index_df["E"].apply(lambda x: x.right)
            index_df.drop(columns=["E"], inplace=True)

            # Convert back to MultiIndex
            df.index = pd.MultiIndex.from_frame(index_df)

        df = df.reset_index()

        # Determine write mode
        mode = "a" if os.path.exists(file) else "w"
        if_sheet_exists = "replace" if mode == "a" else None

        # Write to Excel
        with pd.ExcelWriter(file, mode=mode, engine="openpyxl", if_sheet_exists=if_sheet_exists) as writer:
            df.to_excel(writer, index=False, sheet_name="SMP")

    def iterate_xs_samples(self):
        """
        Iterate samples one by one and shape them as a :obj:`~sandy.xs.Xs`
        dataframe, but with mutligroup structure.
        This output should be passed to :obj:`~sandy.xs.Xs._perturb`.
        The function is called by :obj:`~sandy.endf6..Endf6.apply_perturbations`.

        Yields
        ------
        n : `int`
            sample ID.
        s : `pd.DataFrame`
            dataframe of perturbation coefficients with:
                
                - columns: `pd.MultiIndex` with levels `"MAT"` and `"MT"`
                - index: `pd.IntervalIndex` with multigroup structure

        Notes
        -----
        If samples refer to redundant MT number, the same identical samples
        are passed one level down to the partial MT components.
        For instance:

            - MT=4 samples will be assigned to MT=50-91
            - MT=1 samples will be assigned to MT=2 and MT=3
            - MT=18 samples will be assigned to MT=19-21 and MT=38
        
        .. important:: Assigning samples from redundant MT number to partial
                      components only applies if the partial components do not
                      have their own samples, and it only goes one level deep.

        Examples
        --------

        Get samples fot MT=1.
        >>> import sandy
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10010)
        >>> smps2 = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=2)))[33]

        Copy samples each time to a redundant or partial MT.

        >>> smps1 = sandy.Samples(smps2.data.reset_index().assign(MT=1).set_index(["MAT", "MT", "E"]))
        >>> smps3 = sandy.Samples(smps1.data.reset_index().assign(MT=3).set_index(["MAT", "MT", "E"]))
        >>> smps18 = sandy.Samples(smps1.data.reset_index().assign(MT=18).set_index(["MAT", "MT", "E"]))
        >>> smps19 = sandy.Samples(smps1.data.reset_index().assign(MT=19).set_index(["MAT", "MT", "E"]))
        >>> smps27 = sandy.Samples(smps1.data.reset_index().assign(MT=27).set_index(["MAT", "MT", "E"]))
        >>> smps4 = sandy.Samples(smps1.data.reset_index().assign(MT=4).set_index(["MAT", "MT", "E"]))
        >>> smps51 = sandy.Samples(smps1.data.reset_index().assign(MT=51).set_index(["MAT", "MT", "E"]))
        >>> smps101 = sandy.Samples(smps1.data.reset_index().assign(MT=101).set_index(["MAT", "MT", "E"]))
        >>> smps452 = sandy.Samples(smps1.data.reset_index().assign(MT=452).set_index(["MAT", "MT", "E"]))

        Check that samples are passed correctly to daughter MTs (only one level deep).

        >>> expected = pd.MultiIndex.from_product([[125], [51]], names=["MAT", "MT"])
        >>> assert next(smps51.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [4] + list(sandy.redundant_xs[4])], names=["MAT", "MT"])
        >>> assert next(smps4.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [1] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps1.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [3] + list(sandy.redundant_xs[3])], names=["MAT", "MT"])
        >>> assert next(smps3.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [1] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps1.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [18] + list(sandy.redundant_xs[18])], names=["MAT", "MT"])
        >>> assert next(smps18.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [27] + list(sandy.redundant_xs[27])], names=["MAT", "MT"])
        >>> assert next(smps27.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [101] + list(sandy.redundant_xs[101])], names=["MAT", "MT"])
        >>> assert next(smps101.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [452] + list(sandy.redundant_xs[452])], names=["MAT", "MT"])
        >>> assert next(smps452.iterate_xs_samples())[1].columns.equals(expected)


        In this example the original covariance contains data for MT=1 and MT=51.
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 942400)
        >>> smps = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=[1, 51])))[33]

        Then, since MT=1 is redundant, samples are passed to its partial components (MT=2 and MT=3).
        >>> expected = pd.MultiIndex.from_product([[9440], [1, 51] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps.iterate_xs_samples())[1].columns.equals(expected)
        
        If case one of the partial components already has samples, i.e., MT=2...
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 942400)
        >>> smps = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=[1, 2, 51])))[33]

        Then the MT=1 samples are not passed to the partial components, which 
        in this case it means that MT=2 is not changed and MT=3 is not created.

        >>> expected = pd.MultiIndex.from_product([[9440], [1, 2, 51]], names=["MAT", "MT"])
        >>> assert next(smps.iterate_xs_samples())[1].columns.equals(expected)

        Default use case for MF35.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 942390)
        >>> smps = tape.get_perturbations(2, njoy_kws=dict(err=1, nubar=False, mubar=False))

        Check that output is not empty, and with correct shape.

        >>> assert(next(smps[35].iterate_xs_samples())[1].shape == (240, 5))
        """
        # the tests in this docstrings are very slow...might improve it
        
        levels = Xs._columnsnames
        df = self.data.unstack(level=levels)
        
        # -- Iterate over samples
        for n, p in df.T.groupby(level=self._columnsname):
            s = p.T.droplevel(self._columnsname, axis=1)
            adds = []
            for mat in s.columns.get_level_values("MAT").unique():
                
                # -- Iterate redundant xs (from MT107 to MT1)
                for k, v in redundant_xs.items():
                    if not (mat, k) in s.columns:
                        continue
                    daughters = pd.MultiIndex.from_product([[mat], v], names=["MAT", "MT"])

                    # Only give perturbation for redundant xs to daughters if no perturbation
                    # for partial cross section is found
                    if s.columns.intersection(daughters).empty:
                        
                        # This goes only 1 level deep.
                        # Then, MT=1 perturbations will be given to MT=2 and MT=3
                        # without descending any further
                        add = pd.DataFrame(
                            np.tile(s[(mat, k)].values, (daughters.size, 1)).T,
                            index=s.index,
                            columns=daughters,
                            )
                        adds.append(add)
            s = pd.concat([s, *adds], axis=1)
            yield n, s

    def _std_convergence(self):
        smp = self.data
        rng = range(2, smp.shape[0])
        foo = lambda x: smp.loc[:x].std()
        return pd.DataFrame(map(foo, rng), index=rng)

    def _mean_convergence(self):
        smp = self.data
        rng = range(1, smp.shape[0])
        foo = lambda x: smp.loc[:x].mean()
        return pd.DataFrame(map(foo, rng), index=rng)

    @classmethod
    def from_excel(cls, file, beg=None, end=None):
        """
        Read perturbation coefficients (for nubar and xs) from an Excel file.
        
        The file format should be compatible with what is written in 
        :obj:`~sandy.samples.Samples.to_excel`.
        
        Parameters
        ----------
        file : `str`
            Excel file name.
        beg : `int`
            First sample to consider. Default is first in dataframe.
        end : `int`
            Last samples to conisder. Default is last in dataframe.

        Returns
        -------
        smp : :obj:`~sandy.samples.Samples`
            Samples dataframe.

        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> from sandy.samples import Samples
    
        Create a sample Excel file.

        >>> data = {
        ...     "MT": [102, 102, 18, 18],
        ...     "ELEFT": [1e-5, 1.0, 1e-5, 1.0],
        ...     "ERIGHT": [1.0, 20.0, 1.0, 20.0],
        ...     0: [1.01, 0.99, 1.02, 1.03],
        ...     1: [1.02, 1.00, 1.01, 1.05],
        ... }
        >>> df = pd.DataFrame(data)
        >>> df.to_excel("test_samples.xlsx", sheet_name="SMP", index=False)
    
        Read the data back.

        >>> samples = Samples.from_excel("test_samples.xlsx")
        
        Check the resulting dataframe.

        >>> assert isinstance(samples, Samples)
        >>> assert samples.data.shape == (4, 2)  # 4 rows, 2 samples

        """
        df = pd.read_excel(file, sheet_name="SMP")
        
        # -- Everything before ELEFT is part of the MultiIndex
        loc = df.columns.get_loc("ELEFT")
        idx = df.iloc[:, :loc]

        # -- ELEFT and ERIGHT are combined into an IntervalIndex (multigroup)
        idx.insert(loc, "E", pd.IntervalIndex.from_arrays(df.ELEFT, df.ERIGHT))

        idx = pd.MultiIndex.from_frame(idx)

        df = df.iloc[:, loc+2:].loc[:, beg:end].reset_index(drop=True)
        df.index = idx
        return cls(df)
