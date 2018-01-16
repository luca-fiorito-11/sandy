# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 09:01:48 2017

@author: lfiorito

..From ENDF-6 manual::
    File 35 contains covariance matrices for the energy distribution of 
    secondary particles given in File 5.
    The data in File 5 are normally given in the Laboratory system, and are 
    expressed as normalized probability distributions.

    Each covariance matrix applies to the complete secondary energy 
    distributions for the broad incident energy range specified, regardless 
    of how these secondary energy distributions are specified, or broken down 
    into various components, in File 5.
    No covariances between the different incident energy ranges are allowed.
    Also, no covariances linking different materials or reaction types are 
    allowed.
"""
import logging
import sys
import numpy as np
from sandy.mf33 import Cov33
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard

class MF(MFstandard):
    
    def __init__(self, endf):
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend

    def __setitem__(self, key, item):
        """
        Filter the covariance sections to set in the dictionary according to 
        the following rules:
            * do not set if ``MT`` is not requested
            * do not set if section is empty
        """
        from sandy.sandy_input import options
        if "mt" in options:
            if key not in options['mt']:
                logging.info(item.prefix + "Delete section because MT was not requested")
                return
        if len(item) == 0:
            logging.warn(item.prefix + "Delete section because it is empty")
            return
        super().__setitem__(key, item)
    
    def sorted_items(self):
        r"""
        Generator that yields covariance instances sorted by ``MT`` and 
        covariance subsection.
        
        Outputs:
            - :```key``: :
                (tuple) keyword ``(MT,i)`` where ``i`` is the index of the 
                covariance subsection
            - :``item``: :
                (``mf33.Cov33`` instance) covariance object
        """
        mtlist = sorted(self.keys())
        for mt in mtlist:
            for i,item in enumerate(self[mt]):
                key = (mt,i)
                yield key, item

    def sampling(self):
        r"""
        """
        from sandy.sandy_input import options
        nsmp = options['samples']
        logging.info(self.prefix + "Start sampling procedure")
        if not self:
            logging.error(self.prefix + "No covariance available")
            sys.exit()
        for key,sec in self.sorted_items():
            logging.debug(sec.prefix + "Draw samples for incoming neutron energy {:.3e} <= E <= {:.3e} eV".format(sec.elo, sec.ehi))
            sec.smp = sec.sampling(nsmp)
    
    def write_to_csv(self, filename):
        r"""
        Write covariance matrices to several ``csv`` file, each correspondig to 
        a different ``MT`` reaction and section.
        
        Inputs:
            - :``filename``: :
                (string) name of the ``csv`` file
        """
        from sandy.sandy_input import outdir
        from os.path import join
        namecov = filename + ".cov_mf{}_mt{}_{:.3e}_{:.3e}.csv"
        namesmp = filename + ".smp_mf{}_mt{}_{:.3e}_{:.3e}.csv"
        for (mt,i),sec in self.sorted_items():
            file = join(outdir, namecov.format(self.mf, mt, sec.elo, sec.ehi))
            logging.info("MF{} MT{} EMIN={:.3e} eV EMAX={:.3e} eV : Writing covariance to file  '{}'".format(self.mf, mt, sec.elo, sec.ehi, file))
            sec.cov.write_to_csv(file)
            file = join(outdir, namesmp.format(self.mf, mt, sec.elo, sec.ehi))
            logging.info("MF{} MT{} EMIN={:.3e} eV EMAX={:.3e} eV : Writing samples to file     '{}'".format(self.mf, mt, sec.elo, sec.ehi, file))
            sec.smp.write_to_csv(file)



class MT(list, MTstandard):
    r"""
    ``MT` section for ``MF35``.
    This section contains ``nk`` subsections.
    Each subsection covers a covariance matrix for one incident particle 
    energy range.
    """

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.text = self.read_text(endf, getback=True)
        nk = endf.read_cont([0,1,2,3,5])
        logging.debug(self.prefix + "{} covariance subsections were found".format(nk))
        for k in range(nk):
            section = SubSection(endf)
            self.append(section)

    @property
    def xx(self):
        r"""
        Union energy grid (over all the covariance subsections) for the 
        incident-neutro (:math:`E_{in}`).
        """
        from sandy.records import Grid
        return Grid(np.unique([ sec.domain for sec in self ]))

    def append(self, item):
        r"""
        Append subsection only if the following conditions are verified:
            * the covariance is not empty
        """
        if (item.cov == 0).all():
            logging.info(item.prefix + "Delete subsection because covariance matrix is empty")
            return
        super().append(item)
    
    def get_samples_by_energy(self, energy):
        r"""
        Given the incoming neutron energy return the samples for the 
        corresponding output neutron-energy distribution function.
        
        Inputs:
            - energy :
                (scalar float) incoming neutron energy
        
        Outputs:
            - smp :
                (Samples instance) samples for the output neutron-energy 
                distribution function (raise `error` if samples are not found).
        """
        for section in self:
            if section.contains(energy):
                smp = section.smp
                return smp
        raise AttributeError
    
    def get_cov_by_energy(self, energy):
        r"""
        Given the incoming neutron energy return the covariance for the 
        corresponding output neutron-energy distribution function.
        
        Inputs:
            - energy :
                (scalar float) incoming neutron energy
        
        Outputs:
            - cov :
                (covariance instance) covariance for the output neutron-energy 
                distribution function (raise `error` if covariance is not 
                found).
        """
        for section in self:
            if section.contains(energy):
                cov = section.cov
                return cov
        raise AttributeError
    


class SubSection:
    r"""
    Subsection of a ``MT`` section of ``MF35``.
    This subsection contains the covariance matrix for the energy distributions 
    of secondary particles.
    
    With respect to ``MF33``, a new type of ``lb`` subsection is defined 
    ``(lb=7)``.
    Covariances in ``MF35`` refer to normalized probabilities, therefore it is 
    assumed natural to specify the covariance matrices as absolute covariances 
    of the normalized probabilities rather than the corresponding relative 
    covariances.
    
    The ``lb=7`` subsection is similar to an ``lb=5`` subsection, but with 
    entries that are absolute rather than relative.
    """
    
    def __repr__(self):
        return "<MF{} MT{} CovSec>".format(self.mf, self.mt)

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.elo = endf.line.CONT.c1
        self.ehi = endf.line.CONT.c2
        self.cov = Cov35.read_endf(endf)
    
    @property
    def prefix(self):
        return "MF{} MT{} : ".format(self.mf, self.mt)

    @property
    def domain(self):
        return self.elo, self.ehi
    
    def contains(self, item):
        r"""
        Check docstring of function ``sandy.functions.contains``.
        """
        from sandy.functions import contains
        return contains(item, *self.domain)

    def sampling(self, nsmp):
        from sandy.records import Samples
        smp = self.cov.sampling(nsmp)
        samples = Samples(self.cov.xx, smp)
        return samples



class Cov35(Cov33):

    @staticmethod
    def bin_average(xx, cov):
        r"""
        Normalize covariance matrix dividing by the energy bin.
        
        .. math::
            
            C_{k,k^{'}} = \frac {C_{k,k^{'}}} {\Delta x_{k}  \Delta y_{k^{'}}}
        """
        from sandy.cov import corr2cov
        xx = np.array(xx)
        dx = 1./(xx[1:]-xx[:-1])
        dx = np.insert(dx, len(dx), 0)
        cov = corr2cov(cov, dx)
        return cov
    
    @classmethod
    def read_endf(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in ``lb=7`` 
        format of ``ENDF-6``.
#        ``ENDF-6`` description::
#            Absolute covariance matrix of some cross sections averaged over 
#            some energy intervals.
            
        .. math::
            rcov(x,y) = g(E_x,E_y)

        Inputs:
            - :``endf : 
                (`endf.ENDF` instance)
        """
        from sandy.cov import triu_matrix, Cov
        elo, ehi, ls, lb, nt, ne, data = endf.read_list()
        if lb != 7:
            logging.error(" MF35 : Covariance section must have flag 'lb=7'")
            sys.exit()
        xx = yy = data[:ne]
        f = np.array(data[ne:])
        if ls == 0: # matrix is asymmetric # to be tested
            cov = f.reshape(ne-1, ne-1)
        else: # matrix is symmetric
            cov = triu_matrix(f, ne-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = Cov(np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1))
        cov = cls.bin_average(xx, cov)
        return cls(xx, yy, cov, lb)