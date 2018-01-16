#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:51:47 2017

@author: lfiorito
"""

import sys
import logging
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard


class MF(MFstandard):
    """
    ``MF34`` section of an ``ENDF-6`` file.
    All the corresponding ``MT`` subsections are stored are elements of the 
    dictionary.
    """

    def __init__(self, endf):
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend
    
    @property
    def union_grid(self):
        r"""
        Union grid of the covariance matrices in ``MF34``.
        """
        from sandy.functions import union_grid
        return union_grid(*[ v.xx.tolist() for k,v in self.sorted_items()])

    def __setitem__(self, key, item):
        """
        Filter the covariance sections to set in the dictionary according to 
        the following rules:
            * do not set if ``MT`` is not requested
            * do not set if section is empty
            * do not set if the ``key`` is not in the ``item`` dictionary.
        """
        from sandy.sandy_input import options
        if "mt" in options:
            if key not in options['mt']:
                logging.info(item.prefix + "Delete section because MT was not requested")
                return
        if len(item) == 0:
            logging.warn(item.prefix + "Delete section because it is empty")
            return
        if key not in item:
            logging.warn(item.prefix + "Delete section because covariance for '({0},{0})' was not found".format(key))
            return
        super().__setitem__(key, item)
    
    def sorted_items(self):
        r"""
        Generator that yields covariance instances sorted by ``MF``, ``MT`` 
        and ``MTi``.
        
        Outputs:
            - key :
                (tuple) keyword ``(MT,MTi)``
            - item :
                (``mf33.Cov33`` instance) covariance object
        """
        mtlist = sorted(self.keys())
        for x1 in mtlist:
            mtilist = sorted(self[x1].keys())
            for x2 in mtilist:
                if x2 not in self:
                    logging.info("MF{0} MT{1} MTi{2} : Skip section because covariances for MT{2} were not found".format(self.mf, x1, x2))
                    continue
                else:
                    keys = sorted(self[x1][x2].keys())
                    for (k1,k2) in keys:
                        item = self[x1][x2][k1,k2]
                        key = (x1,k1),(x2,k2)
                        yield key, item

    def sampling(self):
        """
        Extract sample from the union covariance matrix of file `MF34` and 
        sort them by keys.
        Each key is a tuple `(MT,L)`, where `L` is the order of the 
        Legendre Polynomial coefficient.
        """
        from sandy.cov import Ucov
        from sandy.sandy_input import options
        nsmp = options['samples']
        logging.info(self.prefix + "Start sampling procedure")
        data = Ucov({ k : v for k,v in self.sorted_items() })
        if not data:
            logging.error(self.prefix + "No covariance available")
            sys.exit()
        logging.debug(self.prefix + "Draw samples from union covariance matrix")
        for (mt,l),s in data.sampling(nsmp).items():
            self[mt][mt][l,l].smp = s
        if "stdmax" in options:
            self.filter_samples(options['stdmax'])
    
    def filter_samples(self, threshold):
        """
        Delete samples which have standard deviation above a given threshold.
        
        Inputs:
            - :``threshold``: :
                (scalar) threshold of standard deviation
        """
        for ((mt,l),(mti,li)),sec in self.sorted_items():
            if hasattr(sec, "smp"):
                emask, smask = sec.smp.stdmax(threshold)
                logging.info("MF{} MT{} P{} : Stdev > {:.1f}% in {} energy points".format(self.mf, mt, l, threshold*100., len(emask)))
                for x,s in zip(emask, smask):
                    logging.debug(" - E = {:5e} eV --> stdev = {:.1f}%".format(x, s*100.))

    def write_to_csv(self, filename):
        r"""
        Write covariance matrices to several ``csv`` file, each correspondig to 
        a different ``MT`` reaction and a Legendre Polynomial coefficient.
        
        Inputs:
            - :``filename``: :
                (string) name of the ``csv`` file
        """
        from sandy.sandy_input import outdir
        from os.path import join
        namecov = filename + ".cov_mf{}_mt{}_P{}.csv"
        namesmp = filename + ".smp_mf{}_mt{}_P{}.csv"
        for ((mt,l),(mti,li)),sec in self.sorted_items():
            if hasattr(sec, 'smp'):
                file = join(outdir, namecov.format(self.mf, mt, l))
                logging.info("MF{} MT{} P{} : Writing covariance to file  '{}'".format(self.mf, mt, l, file))
                sec.write_to_csv(file)
                file = join(outdir, namesmp.format(self.mf, mt, l))
                logging.info("MF{} MT{} P{} : Writing samples to file     '{}'".format(self.mf, mt, l, file))
                sec.smp.write_to_csv(file)



class MT(dict, MTstandard):
    r"""
    ``MT`` section of ``MF34``.
    This section contains covariances for angular distributions of secondary 
    particles.
    """

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.text = self.read_text(endf, getback=True)
        self.za, self.awr, self.ltt, self.nmt1 = endf.read_cont([2,4])
        for i in range(self.nmt1):
            sec = SubSection(endf)
            self[sec.mt1] = sec

    @property
    def ltt(self):
        r"""
        Flag to specify the representation used, and it may have the following values
        in ``MF34``.
         - :``ltt=0``: :
             all angular distributions are isotropic
         - :``ltt=1`` :
             the data are given as Legendre coefficient covariances as a 
             function of incident energy, starting with :math:`a_1` or 
             higher order coefficients.
         - :``ltt=2``: :
             the data are given as Legendre coefficients covariances as a
             function of incident energy, starting with :math:`a_0`.
             (This information is redundant in the formats, but is
              considered desirable as an alarm flag.)
         - :``ltt=3``: :
             if either :math:`L` or :math`L1=0` anywhere in the section.
        """
        return self._ltt
    
    @ltt.setter
    def ltt(self, ltt):
        if ltt != 1:
            import pytest
            pytest.set_trace()
            logging.error(self.prefix + "SANDY cannot process 'ltt != 1'.")
            sys.exit()
        self._ltt = ltt
    
    @property
    def nmt1(self):
        r"""
        Number of subsections present in ``MF34`` for various 
        :math:`MT1 \\ geq MT`.
        """
        return self._nmt1
    
    @nmt1.setter
    def nmt1(self, nmt1):
        if nmt1 != 1:
            logging.error(self.prefix + "More than one MT reaction ({}) correlated to angular distribution. SANDY cannot process it.")
            sys.exit()
        self._nmt1 = nmt1

    def __setitem__(self, key, item):
        r"""
        Assign dictionary item only if the following conditions are verified:
            * the item does not refer to a cross-material section
            * the item does not refer to a cross-reaction section
            * the item is not empty
        """
        if item.mat1 != 0:
            logging.debug(self.prefix + "SANDY does not process cross-material covariances (MAT{},MATi{})  in section MF34".format(self.mat, item.mat1))
            return
        if item.mt1 != self.mt:
            logging.debug(self.prefix + "SANDY does not process cross-reaction covariances (MT{},MTi{}) in file MF34".format(self.mt, item.mt1))
            return
        if not item:
            logging.debug(self.prefix + "No covariance was found")
            return
        super().__setitem__(key, item)



class SubSection(dict):
    r"""
    Subsection of ``MT`` section of ``MF34``.
    
    The energy-dependent covariance matrices for couples of Legendre 
    coefficients are stored in the dictionary with ``(L,L1)`` tuples as keys, 
    where
     - ``L`` is the Index of the Legendre coefficient for reaction ``MT``
     - ``L1`` is the Index of the Legendre coefficient for reaction ``MT1``
    """
    
    def __repr__(self):
        return "<MF34 subsection MATi={0.mat1} MTi={0.mt1}>".format(self)
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.mat1, self.mt1, self.nl, self.nl1 = endf.read_cont([0,1])
        logging.debug("MF34 MT{} : found covariances with MATi={} MTi={}".format(self.mt, self.mat1, self.mt1))
        for isec in range(self.nss):
            subsec = SubSubSection(endf, self.mat1, self.mt1)
            if self.mat1 == 0 and self.mt1 == self.mt:
                self[subsec.key] = subsec
        self.delete_items()
    
    @property
    def prefix(self):
        return "MF{0.mf} MT{0.mt} MATi{0.mat1} MTi{0.mt1} : ".format(self) 

    @property
    def mat1(self):
        r"""
        ``MAT`` value of the covariance second component.
        """
        return self._mat1

    @mat1.setter
    def mat1(self, mat1):
        _mat1 = 0 if mat1 == self.mat else mat1
        self._mat1 = _mat1

    @property
    def mt1(self):
        r"""
        ``MT`` value of the covariance second component.
        """
        return self._mt1

    @mt1.setter
    def mt1(self, mt1):
        self._mt1 = mt1

    @property
    def nl(self):
        r"""
        Number of Legendre coefficients for which covariance data are given 
        for the reaction ``MT``.
        This value must be the same for each subsection.
        The first coefficient is :math:`a_0` if ``ltt=3``, :math:`a_1 \\geq 1` 
        if ``ltt=1``.
        """
        return self._nl

    @nl.setter
    def nl(self, nl):
        self._nl = nl

    @property
    def nl1(self):
        r"""
        Number of Legendre coefficients for which covariance data are given 
        for reaction ``MT1``.
        """
        return self._nl1

    @nl1.setter
    def nl1(self, nl1):
        self._nl1 = nl1

    @property
    def nss(self):
        r"""
        Number of sub-subsections for a given reaction.
        
        The number of sub-subsections :math:`nss` for a given ``MT1`` is 
        :math:`nl \times nl1`, and they are ordered as
        .. math::
            
            (L,L_i) = (1,1),(1,2),...,(NL,NL_1)
        
        Not all *L*-values need be included.
        When :math:`MT1=MT`, redundancy is avoided by giving each 
        sub-subsection only once, when :math:`L_i \\geq L`.
        In this case
        .. math:
            
            NSS = NL (NL+1)/2.
        """
        if self.mt == self.mt1:
            nss = int(self.nl*(self.nl+1)/2)
        else:
            nss = int(self.nl * self.nl1)
        return nss
    
    def __setitem__(self, key, item):
        r"""
        Set covariance object only if the matrix is not empty.
        """
        if not item.has_cov():
            logging.debug(item.prefix + "Delete section because covariance matrix is empty")
            return
        super().__setitem__(key, item.cov)
    
    def delete_items(self):
        r"""
        Delete dictionary item with key ``(L,L1)`` if item is not present 
        for key ``(L,L)`` or ``(L1,L1)``.
        """
        keys = list(self.keys())
        for (l,l1) in keys:
            if l != l1:
                if (l,l) not in self:
                    logging.info(self.prefix + "Found covariance for (L{0},Li{1}), but not for (L{0},Li{0})".format(l, l1))
                    del self[l,l1]
                    continue
                if (l1,l1) not in self:
                    logging.info(self.prefix + "Found covariance for (L{0},Li{1}), but not for (L{1},Li{1})".format(l, l1))
                    del self[l,l1]
                



class SubSubSection:
    r"""
    Sub-subsection for a given ``MAT`` and ``MT1``.
    """
    
    def __init__(self, endf, mat1, mt1):
        from sandy.mf33 import Cov33
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.mat1 = mat1
        self.mt1 = mt1
        self.l, self.l1, self.lct, self.ni = endf.read_cont([0,1])
        covs = [ Cov33.read_endf(endf) for ini in range(self.ni) ]
        self.covs = [ cov for cov in covs if cov.lb != 0 ]
    
    @property
    def key(self):
        """
        Dictionary key for this sub-subsection.
        """
        return self.l, self.l1
    
    @property
    def prefix(self):
        return "MF{0.mf} MT{0.mt} MATi{0.mat1} MTi{0.mt1} L{0.l} Li{0.l1} : ".format(self) 
        
    @property
    def ni(self):
        r"""
        Number of ``LIST`` records contained in this sub-subsection.
        """
        return self._ni
    
    @ni.setter
    def ni(self, ni):
        if ni != 0:
            logging.debug(self.prefix + "{} 'ni' sections were found".format(ni))
        self._ni = ni
    
    @property
    def covs(self):
        r"""
        List of covariance matrix sections.
        """
        return self._covs
    
    @covs.setter
    def covs(self, covs):
        self._covs = covs
    
    @property
    def cov(self):
        r"""
        Union covariance matrix for this sub-subsection.
        """
        from sandy.mf33 import Cov33
        return Cov33.merge(*self.covs)
    
    def has_cov(self):
        r"""
        Return ``True`` if the covariance matrix exists and is not empty, 
        otherwise ``False``.
        
        Outputs:
            - :``has_cov``:
                (boolean) flag to determine whether section has covariance 
                matrix
        """
        has_cov = False
        if self.covs:
            has_cov = bool(self.cov.any())
        return has_cov