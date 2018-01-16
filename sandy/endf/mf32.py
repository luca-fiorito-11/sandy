#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 18:17:22 2017

@author: lfiorito
"""

import logging
import sys
import numpy as np
from sandy.mf2 import MT as MT2
from sandy.mf2 import BW as BW2
from sandy.mf2 import RM as RM2
from sandy.mf2 import RRR
from sandy.endf import MF as MFstandard
from sandy.cov import Cov


class SrCov(Cov):
    r"""
    """
    
    @property
    def prefix(self):
        return "COV SR : "
    
    @classmethod
    def read_endf(cls, endf):
        r"""
        """
        if endf.line.CONT.n2 == 0: # Breit-Wigner
            obj = cls.read_endf_cont(endf)
        elif endf.line.CONT.n2 == 1: # Reich-Moore
            obj = cls.read_endf_list(endf)
        return obj

    @classmethod
    def read_endf_cont(cls, endf):
        r"""
        """
        dap = endf.read_cont([0,2,3,4,5])
        cov = np.atleast_2d(dap**2)
        return cls(cov)
    
    @classmethod
    def read_endf_list(cls, endf):
        r"""
        """
        _,_,_,_,mls,_,dap = endf.read_list()
        if mls != 1:
            # This was found only in U232 in TENDL-2015. However the format was wrong
            logging.warning("COV SR : 'mls={}', SANDY cannot process uncertainties for l-dependent scattering radius".format(mls))
            cov = np.atleast_2d([0])
        else:
            cov = np.atleast_2d(dap[0]**2)
        return cls(cov)



class BW(RRR):
    r"""
    Energy range in section ``MF32`` where resonance parameters are 
    represented with the Breit-Wigner formalism.
    """
    
    cov_parms = ['er', 'gn', 'gg', 'gf']
    
    def __repr__(self):
        return "<BW RES COV LCOMP{}>".format(self.lcomp)

    @property
    def lcomp(self):
        r"""
        Covariance format flag.
        """
        return self._lcomp
    
    @lcomp.setter
    def lcomp(self, lcomp):
        if lcomp not in [0,1,2]:
            logging.error(self.prefix + "Covariance format 'lcomp={}' is not allowed".format(lcomp))
        logging.debug(self.prefix + "Found covariance matrix in LCOMP{} format".format(lcomp))
        self._lcomp = lcomp

    def __init__(self, endf):
        super().__init__(endf)
        self.isr = self.n2
        if self.isr > 0: # uncertainty provided to the scattering radius
            self.covsr = SrCov.read_endf(endf)
        if self.lcomp == 0:
            covsec = LCOMP0(endf, self._nls)
        elif self.lcomp == 1:
            covsec = LCOMP1(endf, BW)
        elif self.lcomp == 2:
            covsec = LCOMP2(endf, BW)
        if covsec.cov.size != 0:
            self.covsec = covsec
    
    def get_indices(mpar):
        r"""
        """
        indices = [0, 3, 4, 5]
        return indices[:mpar]



class RM(RRR):
    r"""
    Energy range in section ``MF32`` where resonance parameters are 
    represented with the Reich-Moore formalism.
    """
    
    cov_parms = ['er', 'gn', 'gg', 'gfa', 'gfb']

    def __repr__(self):
        return "<RM RES COV LCOMP{}>".format(self.lcomp)

    @property
    def lcomp(self):
        r"""
        Covariance format flag.
        """
        return self._lcomp
    
    @lcomp.setter
    def lcomp(self, lcomp):
        if lcomp not in [1,2]:
            logging.error(self.prefix + "Covariance format 'lcomp={}' is not allowed".format(lcomp))
        logging.debug(self.prefix + "Found covariance matrix in LCOMP{} format".format(lcomp))
        self._lcomp = lcomp
    
    def __init__(self, endf):
        super().__init__(endf)
        self.isr = self.n2
        if self.isr > 0: # uncertainty provided to the scattering radius
            self.covsr = SrCov.read_endf(endf)
        if self.lcomp == 1:
            self.covsec = LCOMP1(endf, RM)
        elif self.lcomp == 2:
            self.covsec = LCOMP2(endf, RM)
    
    def get_indices(mpar):
        r"""
        """
        indices = [0, 2, 3, 4, 5]
        return indices[:mpar]



class LCOMP:
    r"""
    Parent class for ``LCOMP0``, ``LCOMP1`` and ``LCOMP2``.
    """
    
    def __repr__(self):
        return "<MF32 LCOMP{}>".format(self.lcomp)

    @property
    def prefix(self):
        return "MF{} MT{} RRR LCOMP{} : ".format(self.mf, self.mt, self.lcomp)
    
    @property
    def restype(self):
        r"""
        Class that must be used for the resonance parameters.
        """
        from sandy.mf2 import BWres, RMres
        if self.type is BW:
            return BWres
        elif  self.type is RM:
            return RMres
    
    @property
    def parms(self):
        """
        Resonance parameters for which the covariance matrix is given.
        """
        return self.type.cov_parms[:self.mpar]

    @property
    def nrb(self):
        """
        Number of resonances in this block and for which resonance parameter 
        and covariance data are given in this subsection.
        """
        return self._nrb
    
    @nrb.setter
    def nrb(self, nrb):
        logging.debug(self.prefix + "'nrb={}' resonances in this block".format(nrb))
        self._nrb = nrb
    
    @property
    def nsrs(self):
        r"""
        Number of short-range covariance sections.
        """
        return self._nsrs
    
    @nsrs.setter
    def nsrs(self, nsrs):
        if nsrs > 1:
            # found nsrs=3 in n-91-Pa-231.jeff32
            import pytest
            pytest.set_trace()
            logging.error(self.prefix + "SANDY cannot process more than 1 short-range covariance section, found 'nsrs={}'".format(nsrs))
            sys.exit()
        self._nsrs = nsrs

    @property
    def nlrs(self):
        r"""
        Number of long-range covariance sections.
        """
        return self._nlrs
    
    @nlrs.setter
    def nlrs(self, nlrs):
        if nlrs != 0:
            import pytest
            pytest.set_trace()
            logging.error(self.prefix + "SANDY cannot process any long-range covariance section, found 'nlrs={}'".format(nlrs))
            sys.exit()
        self._nlrs = nlrs

    
    @property
    def indices(self):
        r"""
        List of indices of the covariance elements.
        Basically they are the positions of the parameters of a resonance 
        for which covariances are given, repeated into a list as many times as 
        the number of resonances.
        """
        return self.type.get_indices(self.mpar)

    @property
    def list_means(self):
        r"""
        List of mean values corresponding to the covariance elements.
        """
        return np.ravel([means[self.indices] for means in self.means])
    
    @property
    def list_parms(self):
        r"""
        List of parameter names --- e.g. ``er``, ``gn``, ``gg``, ... --- 
        corresponding to the covariance elements.
        """
        return self.nrb * self.parms
    
    @property
    def keys(self):
        r"""
        List of keys corresponding to the covariance elements.
        Each key has the form as in the example below::
            
            RES-1 er=1e5
        """
        keys = []
        ires = 0
        for i, (parm, mean) in enumerate(zip(self.list_parms, self.list_means)):
            if i % self.mpar == 0:
                ires += 1 
            keys.append("RES-{}  {}={}".format(ires, parm, mean))
        return keys

    def sampling(self, nsmp):
        r"""
        Draw samples (centered in zero) from the covariance matrix.
        Then, assign the samples to the corresponding resonance parameter.
        
        Inputs:
            * :``nsmp``: :
                (scalar integer) number of samples to draw
        
        Outputs:
            * :``list_res``: :
                (list) list of resolved resonances with random samples
        """
        list_res = []
        smp = self.cov.sampling(nsmp).reshape(self.nrb, self.mpar, nsmp)
        for parms,samples in zip(self.means, smp):
            rr = self.restype(*parms)
            for key,ss in zip(self.parms, samples):
                rr.__getattribute__(key).smp = ss
                rr.__getattribute__(key).replace_negative()
            list_res.append(rr)
        self.smp = list_res
#        return list_res

    def write_cov_to_csv(self, file):
        r"""
        Write standard deviation and correlation matrix to ``csv`` file.
        
        Inputs: 
            - :``file``: :
                (string) ``csv`` file path + name
        """
        import csv
        from sandy.functions import div0
        todump = [['Resonance parameters', 'Stdev (%)'] + self.keys]
        rstd = np.abs(div0(self.cov.std, self.list_means)*100.)
        for k, std, corr in zip(self.keys, rstd, self.cov.corr*1000.):
            todump.append([k] + [std] + corr.tolist())
        with open(file, 'w') as fout:
            writer = csv.writer(fout, quoting=csv.QUOTE_ALL)
            writer.writerows(todump)

    def write_smp_to_csv(self, file):
        r"""
        Write samples to ``csv`` file.
        
        ..Important::
            This method looks for attribute ``smp`` which is not initialized 
            with the instance.
            Make sure to run method ``sampling`` to initialize attribute 
            ``smp``.
        
        Inputs: 
            - :``file``: :
                (string) ``csv`` file path + name
        """
        import csv
        from sandy.functions import div0
        perts = []
        for i,res in enumerate(self.smp):
            for key in self.parms:
                pert = div0(res.__getattribute__(key).smp, res.__getattribute__(key)) - 1
                perts.append(pert)
        perts = np.array(perts)
        todump = [["Relative perturbations"]]
        todump.append(['Resonance parameters'] + self.keys)
        for i,pert in enumerate(perts.T):
            line = ['smp_{}'.format(i+1)] + (pert*100.).tolist()
            todump.append(line)
        with open(file, 'w') as fout:
            writer = csv.writer(fout, quoting=csv.QUOTE_ALL)
            writer.writerows(todump)



class LCOMP0_L(LCOMP):
    r"""
    Parent  Resolved Resonance Formats (LCOMP=1).
    """

    lcomp = 0
    type = BW

    @property
    def prefix(self):
        return "MF32 MT151 LCOMP0 L{} : ".format(self.l)

    @property
    def means(self):
        """
        Array of mean values for the resonance parameters.
        """
        structured_tab = np.array(self.tab).reshape(self.nrb, 18)
        return structured_tab[:,:6]
    
    @property
    def std(self):
        """
        Array of standard deviations for the resonance parameters.
        """
        return self.cov.std
    
    @property
    def cov(self):
        """
        Covariance matrix.
        There is no correlation for parameters belonging to different 
        resonances.
        """
        from sandy.cov import up2down
        cov = Cov(np.zeros((self.nrb*self.mpar, self.nrb*self.mpar)))
        for i,items in enumerate(np.array(self.tab).reshape(self.nrb, 18)):
            cov[i*self.mpar,i*self.mpar] = items[6]
            cov[i*self.mpar+1,i*self.mpar+1] = items[7]
            cov[i*self.mpar+1,i*self.mpar+2] = items[8]
            cov[i*self.mpar+2,i*self.mpar+2] = items[9]
            cov[i*self.mpar+1,i*self.mpar+3] = items[10]
            cov[i*self.mpar+2,i*self.mpar+3] = items[11]
            cov[i*self.mpar+3,i*self.mpar+3] = items[12]
        return up2down(cov)
    
    @property
    def mpar(self):
        """
        Number of parameters per resonance in this block which have covariance 
        data.
        """
        return 4

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.awri, _, self.l, _, _, self.nrb, self.tab = endf.read_list()
        if self.cov.nnegvar > 0:
            logging.error(self.prefix + "Found {} negative variances".format(self.cov.nnegvar))
            sys.exit()
    


class LCOMP0(list, LCOMP):
    r"""
    Read Compatible Resolved Resonance Format (LCOMP=0).
    """

    lcomp = 0
    type = BW
    
    @property
    def means(self):
        r"""
        Array of mean values for the resonance parameters.
        Loop over the resonance parameters sections for different *l*-values.
        """
        return np.concatenate([ sec.means for sec in self ])
    
    @property
    def std(self):
        r"""
        Array of standard deviations for the resonance parameters.
        Loop over the resonance parameters sections for different *l*-values.
        """
        return np.concatenate([ sec.std for sec in self ])
    
    @property
    def cov(self):
        r"""
        Covariance matrix.
        There is no correlation for parameters belonging to different 
        resonance sections (different *l*-values).
        """
        size = self.std.size
        cov = Cov(np.zeros((size,size)))
        beg = 0
        for sec in self:
            end = beg + sec.std.size
            cov[beg:end, beg:end] = sec.cov
            beg = end
        return cov
    
    @property
    def mpar(self):
        r"""
        Number of parameters per resonance in this block which have covariance 
        data.
        """
        return 4
    
    @property
    def nrb(self):
        r"""
        Number of resonances for which resonance parameters and covariance 
        data are given.
        Sum the resonances in each resonance section.
        """
        return np.sum([ sec.nrb for sec in self ])

    def __init__(self, endf, nls):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        for i in range(nls):
            self.append(LCOMP0_L(endf))



class LCOMP1(LCOMP):
    r"""
    General Resolved Resonance Formats (LCOMP=1).
    """
    
    lcomp = 1
    
    @property
    def cov(self):
        r"""
        Covariance matrix.
        """
        from sandy.cov import triu_matrix
        return Cov(triu_matrix(self.tab[self.nrb*6:], self.mpar*self.nrb))

    @property
    def means(self):
        r"""
        Array of best-estimate values for the resonance parameters.
        """
        return np.array(self.tab[:self.nrb*6]).reshape(self.nrb,6)
    
    @property
    def std(self):
        r"""
        Array of standard deviations for the resonance parameters.
        """
        return self.cov.std

    @property
    def nsrs(self):
        r"""
        Number of subsections with covariances among parameters of 
        specified resonances (short-range).
        """
        return self._nsrs
    
    @nsrs.setter
    def nsrs(self, nsrs):
        if nsrs > 1:
            logging.error(self.prefix + "More than 1 short-range covariance. This option is not yet implemented")
            sys.exit()
        self._nsrs = nsrs

    @property
    def nlrs(self):
        r"""
        Number of subsections containing data on long-range parameter 
        covariances.
        """
        return self._nlrs
    
    @nlrs.setter
    def nlrs(self, nlrs):
        if nlrs > 0:
            logging.error(self.prefix + "Long-range covariances were found. This option is not yet implemented")
            sys.exit()
        self._nlrs = nlrs


    @property
    def mpar(self):
        r"""
        Number of parameters per resonance in this block which have covariance 
        data.
        """
        return self._mpar
    
    @mpar.setter
    def mpar(self, mpar):
        if self.type == BW and mpar == 5:
            logging.error(self.prefix + "SANDY cannot process covariances for the competitive width 'gx'")
            self._mpar = 4
        else:
            self._mpar = mpar
        logging.debug(self.prefix + "'mpar={}' are the parameters per resonance in this block which have covariance data".format(self._mpar))

    @property
    def nvs(self):
        r"""
        Number of covariance elements listed for this block of resonances.
        """
        return int((self.nrb*self.mpar*(self.nrb*self.mpar+1))/2)
    
    def __init__(self, endf, TYPE):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.type = TYPE
        self.awri, self.nsrs, self.nlrs = endf.read_cont([1,2,3])
        _, _, self.mpar, _, _, self.nrb, self.tab = endf.read_list()
        logging.debug(self.prefix + "'nvs={}' covariance elements are listed for this block".format(self.nvs))
        if self.cov.nnegvar > 0:
            logging.error(self.prefix + "Found {} negative variances".format(self.cov.nnegvar))
            sys.exit()




class LCOMP2(LCOMP):
    r"""
    Resolved Resonance Compact Covariance Formats (LCOMP=2).
    """
    
    lcomp = 2
    
    @property
    def nm(self):
        r"""
        Number of lines containing correlations.
        """
        return self._nm
    
    @nm.setter
    def nm(self, nm):
        logging.debug(self.prefix + "'nm={}' lines of correlations records are found".format(nm))
        self._nm = nm
    
    def __init__(self, endf, TYPE):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.type = TYPE
        self.awri, _, _, self.lrx, _, self.nrb, self.tab = endf.read_list()
        self.nnn = endf.line.CONT.l2
        self.nm = endf.line.CONT.n1
        self.corr = Cov(endf.read_compact_cov())

    @property
    def mpar(self):
        r"""
        Number of parameters per resonance in this block which have covariance 
        data.
        """
        from sandy.functions import div0
        mpar = div0(self.nnn, self.nrb)
        if mpar != int(mpar):
            raise NotImplementedError("'nnn' is not divisible by 'nrb'")
        return int(mpar)
    
    @property
    def std(self):
        r"""
        Array of standard deviations for the resonance parameters.
        """
        tab = np.array(self.tab).reshape(self.nrb, 12)
        res = tab[:,6:]
        return np.ravel([ stds[self.indices] for stds in res ])
    
    @property
    def cov(self):
        r"""
        Covariance matrix.
        """
        from sandy.cov import corr2cov
        return corr2cov(self.corr, self.std)
    
    @property
    def means(self):
        r"""
        Array of best-estimate values for the resonance parameters.
        """
        tab = np.array(self.tab).reshape(self.nrb, 12)
        return tab[:,:6]


class MF(MFstandard):
    """
    ``MF32`` section of an ``ENDF-6`` file.
    All the corresponding ``MT`` subsections are stored are elements of the 
    dictionary.
    """
    
    def __init__(self, endf):
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Processing ENDF-6 section MF32 : COVARIANCES OF RESONANCE PARAMETERS "))
        logging.debug(80*"+")
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend
    
#    def __setitem__(self, key, item):
#        """
#        Filter the covariance sections to set in the dictionary according to 
#        the following rules:
#            * do not set if ``MT`` is not requested
#            * do not set if item is empty
#        """
#        from sandy.sandy_input import options
#        if "mt" in options:
#            if key not in options['mt']:
#                logging.info(item.prefix + "Delete section because MT was not requested")
#                return
#        if len(item) == 0:
#            logging.warn(item.prefix + "Delete section because no covariance was found for the RRR")
#            return
#        super().__setitem__(key, item)



class MT(MT2):
    """
    ``MT`` section of file ``MF32``.
    This section is used to store resonance parameters' covariances.
    
    All attributes and methods can be found (with their description) in 
    *mf2.MT*.
    
    .. Note::
        Altough inheriting from ``mf2.MT``, the ``__init__`` method is 
        rewritten.
    """
    
    @property
    def lrf(self):
        r"""
        Copy from ``mf2.MT``.
        """
        return self._lrf

    @lrf.setter
    def lrf(self, lrf):
        if self.lru == 0 and lrf != 0:
            logging.error(self.prefix + "When flag 'lru=0', then 'lrf=0', not 'lrf={}'".format(lrf))
            sys.exit()
        if self.lru == 1:
            if lrf not in [0,1,2,3,4,7]:
                logging.error(self.prefix + "Format 'lrf={}' is not allowed for RRR covariances".format(lrf))
                sys.exit()
            if lrf not in [0,1,2,3]:
                logging.error(self.prefix + "SANDY cannot process RRR covariances in format 'lrf={}'".format(lrf))
                sys.exit()
        if self.lru == 2 and lrf not in [1,2]:
            logging.error(self.prefix + "Format 'lrf={}' is not allowed for URR parameters".format(lrf))
            sys.exit()
        self._lrf = lrf

    RRR_formalisms = { 0 : RRR,
                       1 : BW,
                       2 : BW,
                       3 : RM}
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.text = self.read_text(endf, getback=True)
        super().__init__(endf)

    def append(self, item):
        """
        Append energy sections to the list of a *mf32.MT* according to 
        the following rules:
            * do not append if the energy section does not contain a 
              covariance either for the resonances or for the scattering radius.
        """
        if not hasattr(item, 'covsec') and not hasattr(item, 'covsr'):
            logging.info(item.prefix + "Delete section because no covariance was found")
            return
        super().append(item)

    def write(self, **kwargs):
        r"""
        This is a copy of ``endf.MT.write``.
        """
        from sandy.endf import ENDF
        text = ENDF.add_control(self.text, self.mat, self.mf, self.mt)
        return text

    def write_to_csv(self, name):
        r"""
        Write covariance 
        """
        from sandy.sandy_input import outdir
        from os.path import join
        for i,sec in enumerate(self):
            if not hasattr(sec, "covsec"):
                continue
            file = join(outdir, name + ".cov_mf32_mt151_sec{}_lcomp{}.csv".format(i, sec.lcomp))
            logging.info(self.prefix + "Writing covariance matrices to file '{}'".format(file))
            sec.covsec.write_cov_to_csv(file)
            file = join(outdir, name + ".samples_mf32_mt151_sec{}_lcomp{}.csv".format(i, sec.lcomp))
            logging.info(self.prefix + "Writing samples to file '{}'".format(file))
            sec.covsec.write_smp_to_csv(file)