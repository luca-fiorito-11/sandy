import pdb
import logging

import numpy as np
import pandas as pd

import sandy
from .basecov import BaseCov

__author__ = "Luca Fiorito"
__all__ = [
        "LpcCov",
        ]

class LpcCov(BaseCov):
    """
    Covariance matrix for energy dependent Legendre polynomials coefficients covariance.
    Covariances can be stored for:
        - individual polynomial coefficients,
        - cross polynomial coefficients,
        - cross isotopes,

    Attributes
    ----------
    index : `pandas.MultiIndex`
        index with four levels:
            - `MAT` : `int`, MAT number to identify the isotope
            - `MT` : `int`, MT number to identify the reaction
            - `L` : `int`, polynomial order
            - `E` : `float`, energy of the incident particle
    columns : `pandas.MultiIndex`
        see `index`
    values : `numpy array`
        matrix coefficients
    
    Methods
    -------
    from_endf6
        Extract global cross section/nubar covariance matrix from 
        `sandy.Endf6` instance
    from_list
        Extract global cross section/nubar covariance matrix from iterables 
        of `sandy.formats.utils.EnergyCov` instances
    get_samples
        Extract perturbations from global cross section/nubar covariance matrix
    get_section
        Extract section of global cross section/nubar covariance matrix as a 
        `sandy.formats.utils.EnergyCov` instance
    """
    
    labels = ["MAT", "MT", "L", "E"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = self.labels
        self.columns.names = self.labels
    
    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract global Legendre Polynomials coefficients covariance matrix 
        from `sandy.Endf6` instance.
        
        Parameters
        ----------
        endf6 : `sandy.Endf6`
            ENDF-6 file instance containing covariance sections
        
        Returns
        -------
        `sandy.LpcCov`
            covariance matrix for Legendre polynomial coefficients found in ENDF-6 file instance
        """
        tape = endf6.filter_by(listmf=[34])
        data = []
        # Loop MF/MT
        logging.debug("found {} covariance sections".format(len(tape)))
        for (mat,mf,mt), text in tape.TEXT.iteritems():
            X = tape.read_section(mat, mf, mt)
            # Loop subsections
            logging.debug("reading section MAT={}/MF={}/MT={}".format(mat, mf, mt))
            logging.debug("found {} subsections".format(len(X["REAC"])))
            for (mat1,mt1), rsec in X["REAC"].items():
                if mat1 == 0:
                    mat1 = mat
                logging.debug("\treading subsection MAT1={}/MT1={}".format(mat1, mt1))
                logging.debug("\tfound {} P sub-subsection".format(len(rsec["P"])))
                for (l,l1), psec in rsec["P"].items():
                    logging.debug("\treading sub-subsection for (P{},P{})".format(l,l1))
                    logging.debug("\tfound {} NI-type sub-sub-subsection".format(len(psec["NI"])))
                    covs = []
                    for i,nisec in psec["NI"].items():
                        logging.debug("\t\treconstruct covariance from NI-type section LB={}".format(nisec["LB"]))
                        if nisec["LB"] == 5:
                            foo = sandy.EnergyCov.from_lb5_asym if nisec["LS"] == 0 else sandy.EnergyCov.from_lb5_sym
                            cov = foo(nisec["EK"], nisec["FKK"])
                            covs.append(cov)
                        elif nisec["LB"] == 1:
                            cov = sandy.EnergyCov.from_lb1(nisec["EK"], nisec["FK"])
                            covs.append(cov)
                        elif nisec["LB"] == 2:
                            cov = sandy.EnergyCov.from_lb2(nisec["EK"], nisec["FK"])
                            covs.append(cov)
                        elif nisec["LB"] == 6:
                            cov = sandy.EnergyCov.from_lb6(nisec["EK"], nisec["EL"], nisec["FKL"])
                            covs.append(cov)
                        else:
                            logging.warn("skip LB={} covariance for [({}/{}), ({}/{})]".format(nisec["LB"], mat, mt, mat1, mt1))
                            continue
                    if len(covs) == 0:
                        logging.debug("\tsubsection MAT1={}/MT1={} did not provide accetable covariances".format(mat1, mt1))
                        continue
                    cov = sandy.EnergyCov.sum_covs(*covs)
                    if cov.all().all():
                        logging.warn("\tempty covariance for [({}/{}), ({}/{})]".format(mat, mt, mat1, mt1))
                        continue
                    data.append([(mat, mt, l), (mat1, mt1, l1), cov])
        if not data:
            logging.warn("no lpc covariance was found")
            return pd.DataFrame()
        return cls._from_list(data)
    
    def plot_std(self, display=True, **kwargs):
        """
        Plot standard deviations with seaborn.
        
        Parameters
        ----------
        display : `bool`
            flag to display figure to screen
        
        kwargs : keyword arguments
            extra arguments to pass to `seaborn.lineplot`
        
        Returns
        -------
        `matplotlib.pyplot.Axes`
        """
        std = self.get_std()*100
        df = std.to_frame().reset_index()
        df["L"] = df["L"].astype("category")
        palette = list(colors.keys())[:len(df.L.unique())]
        ax = sns.lineplot(data=df, drawstyle="steps-post", x="E", y="STD", hue="L", palette=palette, style="MT", **kwargs)
        ax.set_xscale("log")
        if (df.STD > 200).any():
            ax.set_yscale("log")
        ax.set(xlabel='energy (eV)', ylabel='stdev (%)')
        if display:
            plt.grid()
            plt.show()
            plt.close()
        return ax

    def filter_p(self, p):
        """
        Delete covariances for Legendre polynomial coefficients with order higher than `p`.
        
        Parameters
        ----------
        p : `int`
            maximum order of Legendre polynomial coefficients
        
        Returns
        -------
        `sandy.LpcCov`
            reduced covariance matrix with Legendre polynomial coefficients of order lower than `p` 
        """
        mask = self.index.get_level_values("L") <= p
        lpccov = self.iloc[mask, mask]
        return LpcCov(lpccov)

    def get_samples(self, nsmp, **kwargs):
        """
        Draw samples from probability distribution centered in 1 and with
        relative covariance in `LpcCov` instance.
        
        Parameters
        ----------
        nsmp : `int`
            number of samples
        
        Returns
        -------
        `sandy.LpcSamples`
        """
        smp = self.to_matrix().sampling(nsmp) + 1
        frame = pd.DataFrame(smp, index=self.index, columns=range(1,nsmp+1))
        if "eig" in kwargs:
            if kwargs["eig"] > 0:
                eigs = cov.eig()[0]
                idxs = np.abs(eigs).argsort()[::-1]
                dim = min(len(eigs), kwargs["eig"])
                eigs_smp = Cov(np.cov(frame.values)).eig()[0]
                idxs_smp = np.abs(eigs_smp).argsort()[::-1]
                print("MF34 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
                diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
                E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
                print("\n".join(E))
        return frame