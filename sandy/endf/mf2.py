# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 15:42:21 2017

@author: lfiorito
"""
import sys
import logging
import numpy as np
from sandy.endf import MF as MFstandard
from sandy.endf import MT as MTstandard
from sandy.records import Tab1
    

class ScalarFloatSmp(float):
    """
    General class for scalar float that are supposed to have samples.
    This class inherits from a float scalar and adds attribute ``smp`` to store 
    random samples centered on the float value.
    """
    
    @property
    def smp(self):
        r"""
        Absolute samples centered on the ``self`` itself.
        """
        return self._smp
    
    @smp.setter
    def smp(self, smp):
        if np.ndim(smp) != 1:
            raise NotImplementedError("Samples of a scalar value must have 'ndim=1', not 'ndim={}'".format(np.ndim(smp)))
        _smp = np.ravel(smp).astype(float)
        if self == 0:
            self._smp = np.zeros_like(_smp) + self
        else:
            self._smp = _smp + self
    
    @property
    def mean(self):
        r"""
        Mean of the scalar samples.
        """
        if not hasattr(self, '_smp'):
            raise NotImplementedError("Attribute 'smp' of 'ScalarFloatSmp' instance has not yet been allocated")
        return self.smp.mean()

    @property
    def std(self):
        r"""
        Standard deviation of the scalar samples.
        """
        if not hasattr(self, '_smp'):
            raise NotImplementedError("Attribute 'smp' of 'ScalarFloatSmp' instance has not yet been allocated")
        return self.smp.std()

    def __getitem__(self, ismp):
        r"""
        Get sample if present, otherwise return best estimate value.
        
        Inputs:
            - :``ismp``: :
                (scalar integer) array-index of the requested sample
        
        Outputs:
            - :``value``: : 
                (``ScalarFloatSmp`` instance) requested sample
        """
        condition = ismp is not None and hasattr(self, '_smp')
        sample = self.smp[ismp] if condition else self
        value = self.__class__(sample)
        return value
    
    def replace_negative(self):
        r"""
        Replace negative samples with best-estimate value.
        """
        sign = np.sign(self)
        mask = self.smp > 0 if sign == -1 else self.smp < 0
        self.smp[mask] = self



class ResRes:
    r"""
    Class with common properties and methods of resolved resonances.
    """
    
    @property
    def er(self):
        return self._er
    
    @er.setter
    def er(self, er):
        self._er = ScalarFloatSmp(er)

    @property
    def aj(self):
        return self._aj
    
    @aj.setter
    def aj(self, aj):
        self._aj = ScalarFloatSmp(aj)

    @property
    def gn(self):
        return self._gn
    
    @gn.setter
    def gn(self, gn):
        self._gn = ScalarFloatSmp(gn)

    @property
    def gg(self):
        return self._gg
    
    @gg.setter
    def gg(self, gg):
        self._gg = ScalarFloatSmp(gg)

    @property
    def gf(self):
        return self._gf
    
    @gf.setter
    def gf(self, gf):
        self._gf = ScalarFloatSmp(gf)

    @property
    def gfa(self):
        return self._gfa
    
    @gfa.setter
    def gfa(self, gfa):
        self._gfa = ScalarFloatSmp(gfa)

    @property
    def gfb(self):
        return self._gfb
    
    @gfb.setter
    def gfb(self, gfb):
        self._gfb = ScalarFloatSmp(gfb)

    def is_same(self, res):
        r"""
        Given a resolved resonance of type ``mf2.ResRes``, check if it 
        corresponds to the same resonance as ``self``.
        
        The check is successful if they have the same spin and resonance 
        energy, the latter within a tolerance of :math:`10^{-5}%`.
        
        Inputs:
            * :``res``: :
                (``mf2.ResRes`` instance) the resolved resonance to compare.
        
        Outputs:
            * :``same``: :
                (boolean) ``True`` if they are the same resonance, else 
                ``False``
        """
        if res.parms != self.parms:
            raise NotImplementedError("To compare resonances they must have the same parameters 'parms'")
        same = np.isclose(res.er, self.er) and res.aj == self.aj
        if same:
            for parm in self.parms:
                if not np.isclose(self.__getattribute__(parm), res.__getattribute__(parm)):
                    logging.warn(self.prefix + "'{0}={1:.6e}' in MF2, '{0}={2:.5e}' in MF32".format(parm, self.__getattribute__(parm), res.__getattribute__(parm)))
        return same

    def tolist(self, ismp=None):
        r"""
        Return the resolved resonance parameters into a list.
        If the index of a sample is specified, return the corresponding sample.
        
        Inputs:
            - :``ismp``: :
                (scalar integer) index of the sample
            
        Outputs:
            - :``LIST``: :
                (list) list of parameters
        """
        LIST = [ self.__getattribute__(parm)[ismp] for parm in self.parms ]
        return LIST



class BWres(ResRes):
    r"""
    Breit-Wigner resonance.
    """

    def __repr__(self):
        return "<BW resonance>"
    
    @property
    def prefix(self):
        return "BW RES E={:.3e} SPI={:.1f} : ".format(self.er, self.aj)
    
    parms = 'er', 'aj', 'gt', 'gn', 'gg', 'gf'

    @property
    def gt(self):
        r"""
        Resonance total width :math:`\Gamma` as the sum of the 
        partial widths

        .. math::
            \Gamma = \Gamma_{n} + \Gamma_{g} + \Gamma_{f} + \Gamma_{x}
        """
        smp = 0
        for g in self.gn, self.gg, self.gf, self.gx:
            smp += g.smp if hasattr(g, 'smp') else g
        if np.ndim(smp) == 1:
            self._gt.smp = smp - self._gt
        return self._gt
    
    @gt.setter
    def gt(self, gt):
        SUM = self.gn + self.gg + self.gf
        if not np.isclose(gt, SUM) and self.lrx == 0:
            logging.warn(self.prefix + "Total width not equal to sum of partial widths: {} != {}".format(gt, SUM))
        self._gt = ScalarFloatSmp(gt)

    def __init__(self, er, aj, gt, gn, gg, gf, lrx=0):
        self.er = er
        self.aj = aj
        self.gn = gn
        self.gg = gg
        self.gf = gf
        self.lrx = lrx
        self.gt = gt
        self.gx = ScalarFloatSmp(gt - gn - gg - gf) if self.lrx == 1 else ScalarFloatSmp(0) # competitive width
        
    def _get_smp(self, ismp):
        r"""
        """
        parms = []
        for parm in self.parms:
            attr = self.__getattribute__(parm)
            if ismp is not None and hasattr(attr, 'smp'):
                parms.append(attr.smp[ismp])
            else:
                parms.append(attr)
        # In case of competitive width
        if self.lrx == 1:
            if ismp is not None and hasattr(self.gx, 'smp'):
                gx = self.gx.smp[ismp]
            else:
                gx = self.gx
            parms[2] = sum(parms[3:], gx)
        return self.__class__(*parms, self.lrx)



class RMres(ResRes):
    r"""
    Reich-Moore resonance.
    """

    def __repr__(self):
        return "<RM resonance>"
    
    @property
    def prefix(self):
        return "RM RES E={:.6e} SPI={:.1f} : ".format(self.er, self.aj)

    parms = 'er', 'aj', 'gn', 'gg', 'gfa', 'gfb'

    def __init__(self, er, aj, gn, gg, gfa, gfb):
        self.er = er
        self.aj = aj
        self.gn = gn
        self.gg = gg
        self.gfa = gfa
        self.gfb = gfb

    def _get_smp(self, ismp):
        r"""
        """
        parms = []
        for parm in self.parms:
            attr = self.__getattribute__(parm)
            if ismp is not None and hasattr(attr, 'smp'):
                parms.append(attr.smp[ismp])
            else:
                parms.append(attr)
        return self.__class__(*parms)



class SR(Tab1):
    r"""
    Class for the energy-dependent scattering radius.
    It inherits from *records.Tab1* since the scattering radius is simply a
    1d tabulated function.
    """
    
    def __repr__(self):
        return "<SR>"



class URR(list):
    """
    Unresolved resonance region (URR) sub-section of ``MF2``.
    
    A ``URR`` instance contains the unresolved resonance sections as a list of 
    strings.
    """
        
    @classmethod
    def read_endf(cls, endf):
        """
        Read the rest of ``ENDF-6`` file.
        """
        text = []
        while True:
            mat, mf, mt, ns = endf.read_control()
            if mf not in [2,32] or mt != 151:
                break
            line = endf.next()
            text.append(line.HEAD)
        return cls(text)

    def write(self, ismp=None):
        r"""
        Return the URR section as a list of strings.
        
        Inputs:
            - ismp :
                (integer) number of the sample

        Ouputs:
            - :``self``: :
                list of strings (the instance is already a list)
        """
        return self



class RRR(list):
    """
    Resonance section for a given energy range.
    
    Found in:
        - Fm255 (jeff32)
    """
    
    def __repr__(self):
        return "<Only Scattering Radius>".format(self.el, self.eh)

    @property
    def prefix(self):
        return "MF{} MT{} RRR : ".format(self.mf, self.mt)

    @property
    def covsec(self):
        r"""
        Resonance parameter covariance matrix.
        """
        return self._covsec
    
    @covsec.setter
    def covsec(self, covsec):
        r"""
        If covariance is empty, do not allocate.
        """
#        from sandy.cov import Cov
#        import pytest
#        pytest.set_trace()
#        if not isinstance(covsec, Cov):
#            raise NotImplementedError("'covsec' must be of type {}, not {}".format(Cov, type(covsec)))
        self._covsec = covsec

    @property
    def covsr(self):
        r"""
        Scattering radius covariance matrix.
        Since we do not process energy-dependent covariances for the 
        scattering radii, it is always a 1 x 1 matrix.
        """
        return self._covsr
    
    @covsr.setter
    def covsr(self, covsr):
        r"""
        If covariance is empty or zero, do not allocate.
        """
        from sandy.cov import Cov
        if not isinstance(covsr, Cov):
            raise NotImplementedError("'covsr' must be of type {}, not {}".format(Cov, type(covsr)))
        if (covsr == 0).all():
            logging.debug(covsr.prefix + "All elements of covariance matrix are zeros")
        else:
            self._covsr = covsr

    @property
    def el(self):
        r"""
        Lower limit for an energy range.
        """
        return self._el

    @el.setter
    def el(self, el):
        self._el = el

    @property
    def eh(self):
        r"""
        Upper limit for an energy range.
        """
        return self._eh

    @eh.setter
    def eh(self, eh):
        r"""
        Check that upper limit for an energy range is larger than the lower 
        limit.
        """
        logging.debug(self.prefix + "Found resonance region defined for {:.3e} eV <= E <= {:.3e} eV".format(self.el, eh) )
        if eh < self.el:
            logging.error(self.prefix + "Lower energy limit is larger than upper energy limit: {:.3} > {:.3}".format(self.el, eh))
            sys.exit()
        elif eh == self.el:
            logging.warn(self.prefix + "Upper energy limit coincides with lower energy limit: {:.3} = {:.3}".format(self.el, eh) )
            sys.exit()
        self._eh = eh
    
    @property
    def lru(self):
        r"""
        Same as in ``mf2.MT``.
        
        .. Note::
            There is no need to check ``lru`` here.
            The check is done in ``mf2.MT``.
            If you are here it is because ``mf2.MT`` sent you.
        """
        return self._lru

    @lru.setter
    def lru(self, lru):
        self._lru = lru

    @property
    def lrf(self):
        r"""
        .. Important::
            All the checks for compatibility between ``lru`` and ``lrf`` are 
            already done in *mf2.MT*.
        """
        return self._lrf

    @lrf.setter
    def lrf(self, lrf):
        if self.lru == 0:
            if lrf == 0:
                logging.debug(self.prefix + "Found scattering radius section 'lrf={}'".format(lrf))
        if self.lru == 1:
            if lrf == 1:
                logging.debug(self.prefix + "Found Single-level Breit Wigner section 'lrf={}'".format(lrf))
            elif lrf == 2:
                logging.debug(self.prefix + "Found Multi-level Breit Wigner section 'lrf={}'".format(lrf))
            elif lrf == 3:
                logging.debug(self.prefix + "Found Reich-Moore section 'lrf={}'".format(lrf))
        self._lrf = lrf
    
    @property
    def nro(self):
        r"""
        Flag designating possible energy dependence of the scattering radius.
        
        - :``nro=0``: :
            radius is energy independent
        - :``nro=1``: :
            radius expressed as a table of energy, radius pairs 
            (e.g. Tm168 endf/b-vii.1)
        """
        return self._nro

    @nro.setter
    def nro(self, nro):
        if nro == 1:
            logging.debug(self.prefix + "Found energy-dependent scattering radius")
        elif nro != 0:
            logging.error(self.prefix + "Flag 'nro' must be either 0 or 1, not {}".format(nro))
            sys.exit()
        if self.lru == 0 and nro != 0:
            logging.error(self.prefix + "For a section where only the scattering radius is given, flag 'nro' must be 0, not 'nro={}'".format(nro))
            sys.exit()
        if nro != 0 and self.mf == 32:
            logging.error(self.prefix + "SANDY cannot process covariances for energy dependent scattering radius")
            import pytest
            pytest.set_trace()
            sys.exit()
        self._nro = nro
    
    @property
    def naps(self):
        r"""
        Flag controlling the use of the two radii, the channel radius ``a`` and 
        the scattering radius ``ap``.
        """
        return self._naps

    @naps.setter
    def naps(self, naps):
        if self.lru == 0 and naps != 0:
            logging.error(self.prefix + "For a section where only the scattering radius is given, flag 'naps' must be 0, not 'naps={}'".format(naps))
            sys.exit()
        self._naps = naps
    
    @property
    def ap(self):
        r"""
        Scattering radius in units of :math:`10^{-12}` cm.
        For ``lrf=1`` to ``lrf=4``, it is assumed to be independent of the 
        channel quantum numbers.
        """
        return self._ap
    
    @ap.setter
    def ap(self, ap):
        logging.debug(self.prefix + "Scattering radius is 'ap={}'".format(ap))
        self._ap = ScalarFloatSmp(ap)

    @property
    def nls(self):
        r"""
        Number of *l*-values (neutron orbital angular momentum) in this 
        energy region.
        A set of resonance parameters is given for each *l*-value.
        """
        return len(self)
    
    @nls.setter
    def nls(self, nls):
        nlsmax = 4
        if nls > nlsmax:
            logging.warn(self.prefix + "The number of l-values 'nls={}' exceeds the maximum permitted 'nlsmax={}'".format(nls, nlsmax))
        if self.lru == 0 and nls != 0:
            logging.error(self.prefix + "For a section where only the scattering radius is given, the number of l-values must be 0, not 'nls={}'".format(nls))
            sys.exit()
        logging.debug(self.prefix + "'nls={}' l-values found in this energy region".format(nls))
        self._nls = nls

    @property
    def nlsc(self):
        """
        Number of *l*-values which must be used to converge the calculation 
        with respect to the incident *l*-value in order to obtain accurate 
        elastic angular distributions.
        """
        return self._nlsc

    @nlsc.setter
    def nlsc(self, nlsc):
        if self._nls > nlsc:
            logging.warn(self.prefix + "Found 'nls > nlsc': {} > {}".format(self.nls, nlsc))
        self._nlsc = nlsc

    @property
    def isr(self):
        """
        Flag to to indicate the presence or absence of the uncertainty data 
        for the scattering radius, meaning:
            - :``isr=0``: :
                no scattering radius uncertainty data
            - :``isr=1``: :
                scattering radius uncertainty data are present
        """
        return self._isr
    
    @isr.setter
    def isr(self, isr):
        if isr == 1:
            logging.debug(self.prefix + "Scattering radius uncertainty data are present, 'isr={}'".format(isr))
        elif isr != 0:
            logging.error(self.prefix + "Flag 'isr' must be either 0 or 1, not 'isr={}'".format(isr))
        self._isr = isr

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.el, self.eh, self.lru, self.lrf, self.nro, self.naps = endf.read_cont()
        if self.nro != 0:
            self.sr = SR.read_endf(endf)
        self.spi, self.ap, self.lad, self.lcomp, self.nls, self.n2 = endf.read_cont()

    def write(self, ismp=None):
        r"""
        Write the instance data into a list of strings according to the 
        ``ENDF-6`` format rules.
        
        Ouputs:
            - :``text``: :
                list of strings
        """
        from sandy.endf import ENDF
        text = []
        text += ENDF.write_cont(self.el, self.eh, self.lru, self.lrf, self.nro, self.naps)
        if self.nro != 0:
            text += self.sr.write_endf()
        text += ENDF.write_cont(self.spi, self.ap[ismp], self.lad, self.lcomp, self.nls, self.n2)
        return text

    def replace_resonance(self, res):
        r"""
        Given an individual resonance section, check if this resonance exists 
        amongst the instance subsections.
        If it does, then replace the existing section with the given one.
        
        Inputs:
            - :``res``: :
                (``RMres`` or ``BWres`` instance) resonance section
        
        Outputs:
            - :``found``: :
                (boolean) ``True`` if the replacement succeded, otherwise 
                ``False``
        """
        if not isinstance(res, self.restype):
            raise NotImplementedError("Input resonance is in wrong format")
        found = False
        for lsec in self:
            found = lsec.replace_resonance(res)
            if found:
                return found
        return found



class BW(RRR):
    r"""
    Energy range where resonance parameters are represented with the 
    Breit-Wigner formalism.
    
    Found in:
        - Bi209 jeff32
    """
    restype = BWres
    
    def __repr__(self):
        return "<BW section>".format(self)
    
    def __init__(self, endf):
        super().__init__(endf)
        for l in range(self._nls):
            sec = L_section_BW(endf)
            self.append(sec)

    def write(self, ismp=None):
        r"""
        Write the instance data into a list of strings according to the 
        ``ENDF-6`` format rules.
        
        Inputs:
            - ismp :
                (integer) number of the sample

        Ouputs:
            - :``text``: :
                list of strings
        """
        text = super().write(ismp=ismp)
        for lsec in self:
            text += lsec.write(ismp=ismp)
        return text



class RM(RRR):
    r"""
    Energy range where resonance parameters are represented with the 
    Reich-Moore formalism.
    """
    restype = RMres

    def __repr__(self):
        return "<RM section>".format(self)

    def __init__(self, endf):
        super().__init__(endf)
        self.nlsc = self.n2
        for l in range(self._nls):
            sec = L_section_RM(endf)
            logging.debug(sec.prefix + "Scattering radius is 'apl={}'".format(sec.apl))
            self.append(sec)

    def write(self, ismp=None):
        r"""
        Write the instance data into a list of strings according to the 
        ``ENDF-6`` format rules.
        
        Inputs:
            - ismp :
                (integer) number of the sample

        Outputs:
            - :``text``: :
                list of strings
        """
        text = super().write(ismp=ismp)
        for lsec in self:
            text += lsec.write(ismp=ismp)
        return text



class L_section(list):
    r"""
    Class of common properties and methods to ``mf2.L_section_BW`` and 
    ``mf2.L_section_RM``.
    """
    
    @property
    def nrs(self):
        r"""
        Number of resolved resonances for a given *l*-value.
        """
        return len(self)

    @nrs.setter
    def nrs(self, nrs):
        r"""
        Check if number of resolved resonances exceeds maximum value allowed.
        """
        nrsmax = 5000
        if nrs > nrsmax:
            logging.warn(self.prefix + "The number of resolved resonances 'nrs={}' for the given l-value 'l={}' exceeds the maximum permitted 'nrsmax={}'".format(nrs, self.l, nrsmax))
        logging.debug(self.prefix + "'nrs={}' resolved resonances are given for this l-value".format(nrs))
        self._nrs = nrs    
    
    def replace_resonance(self, res):
        r"""
        Given a resonance section, check if this resonance exists amongst 
        the instance subsections.
        If it does exist, then replace the existing section with the given one.
        
        Inputs:
            - :``res``: :
                (``RMres`` or ``BWres`` instance) resonance section
        
        Outputs:
            - :``found``: :
                (boolean) ``True`` if the replacement succeded, otherwise 
                ``False``
        """
        found = False
        for i in range(len(self)):
            found = self[i].is_same(res)
            if found:
                self[i] = res
                return found
        return found



class L_section_BW(L_section):
    r"""
    Subsection of BW section for a specific orbital angular momentum *l*.
    """
    
    def __repr__(self):
        return "<BW section for l={0.l}>".format(self)
    
    @property
    def prefix(self):
        return "MF{} MT{} BW L={} : ".format(self.mf, self.mt, self.l)
    
    @property
    def lrx(self):
        r"""
        Flag indicating whether this energy range contains a competitive width:
            * :``lrx=0``: :
                no competitive width is given
            * :``lrx=1``: :
                a competitive width is given, and is an inelastic process to 
                the first excited state.
        """
        return self._lrx
    
    @lrx.setter
    def lrx(self, lrx):
        if lrx == 1:
            logging.debug(self.prefix + "This range contains a competitive width")
        elif lrx != 0:
            logging.error(self.prefix + "Flag 'lrx' must be either 0 or 1, not 'lrx={}'".format(lrx))
        self._lrx = lrx
    
    @property
    def qx(self):
        r"""
        *Q*-value to be added to the incident particle's center-of-mass energy 
        to determine the channel energy for use in the penetrability factor.
        """
        if self.lrx == 0:
            return 0.0
        else:
            return self._qx
    
    @qx.setter
    def qx(self, qx):
        self._qx = qx
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.awri, qx, self.l, self.lrx, _, self.nrs, tab = endf.read_list()
        self.qx = qx
        for i in range(self._nrs):
            res = BWres(*tab[i*6:i*6+6])
            self.append(res)

    def write(self, ismp=None):
        r"""
        Write list of resonance parameters for the given l-value using 
        the BW formalism expressed in the `ENDF-6` format.
        
        Inputs:
            - ismp :
                (integer) number of the sample (standard is not to plot samples, 
                but best estimate)
        
        Outputs:
            - text :
                (list of strings)
        """
        from sandy.endf import ENDF
        tab = np.ravel([ res.tolist(ismp) for res in self ])
        text = ENDF.write_list(self.awri, self.qx, self.l, self.lrx, self.nrs, tab)
        return text



class L_section_RM(L_section):
    r"""
    Subsection of RM section for a specific orbital angular momentum *l*.
    """

    def __repr__(self):
        return "<RM section for l={0.l}>".format(self)

    @property
    def prefix(self):
        return "MF{} MT{} RM L{} : ".format(self.mf, self.mt, self.l)

    @property
    def apl(self):
        """
        *l*-dependent scattering radius.
        If zero, use ``apl=ap``.
        """
        return self._apl
    
    @apl.setter
    def apl(self, apl):
        self._apl = ScalarFloatSmp(apl)

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.awri, self.apl, self.l, _, _, self.nrs, tab = endf.read_list()
        for i in range(self._nrs):
            res = RMres(*tab[i*6:i*6+6])
            self.append(res)

    def write(self, ismp=None):
        r"""
        Write list of resonance parameters for the given l-value using 
        the RM formalism expressed in the ``ENDF-6`` format.
        
        Inputs:
            - ismp :
                (integer) number of the sample (standard is not to plot 
                samples, but best estimate)
        
        Outputs:
            - text :
                (list of strings)
        """
        from sandy.endf import ENDF
        tab = np.ravel([ res.tolist(ismp) for res in self ])
        text = ENDF.write_list(self.awri, self.apl[ismp], self.l, 0, self.nrs, tab)
        return text



class RML(list):
    r"""
    Energy range where resonance parameters are represented with the 
    R-matrix limited formalism.
    
    Found in:
        - Cu63 (jeff32)
    """    

    def __repr__(self):
        return "<RML section>"
    
    def __init__(self, endf):
        self.el, self.eh, self.lru, self.lrf, self.nro, self.naps = endf.read_cont()
        self.ifg, self.krm, njs, self.krl = endf.read_cont([0,1])
        _, _, npp, _, _, _, tab = endf.read_list()
        self.pp = [ ParticlePair(*tab[i*12:12+i*12]) for i in range(npp) ]
        for i in range(njs):
            sec = J_group(endf)
            self.append(sec)
    
    @property
    def krm(self):
        r"""
        Flag to specify which formulae for the R-matrix are to be used.
        - krm=1 :
            single-level Breit-Wigner
        - krm=2 :
            multilevel Breit-Wigner
        - krm=3 :
            Reich Moore
        - krm=4 :
            full R-matrix
        """
        return self._krm

    @krm.setter
    def krm(self, krm):
        if krm != 3:
            logging.error("MF2 MT151 : SANDY cannot handle R-matrix limited format with 'krm={}'".format(krm))
            sys.exit()
        self._krm = krm

    @property
    def ifg(self):
        r"""
        - Ifg=0 :
            `gam` is the channel width in eV
         -ifg=1 :
             `gam` is the reduced-width amplitude in eV:math:^{1/2}:
        """
        return self._ifg

    @ifg.setter
    def ifg(self, ifg):
        self._ifg = ifg
    
    @property
    def krl(self):
        r"""
        - krl=0 :
            non-relativistic kinematics
        - krl=1 :
            relativistic kinematics
        """
        return self._krl

    @krl.setter
    def krl(self, krl):
        self._krl = krl
    
    def write(self, ismp=None):
        from sandy.endf import ENDF
        text = []
        text += ENDF.write_cont(self.el, self.eh, self.lru, self.lrf, self.nro, self.naps)
        text += ENDF.write_cont(0, 0, self.ifg, self.krm, len(self), self.krl)
        tab = np.ravel([ pp.tab for pp in self.pp ])
        text += ENDF.write_list(0, 0, len(self.pp), 0, len(self.pp)*2, tab)
        for sec in self:
            text += sec.write(ismp=ismp)
        return text


class ParticlePair:
    
    def __repr__(self):
        return "<MT{}>".format(self.mt)

    def __init__(self, ma, mb, za, zb, ia, ib, q, pnt, shf, mt, pa, pb):
        self.ma = ma
        self.mb = mb
        self.za = za
        self.zb = zb
        self.ia = ia
        self.ib = ib
        self.q = q
        self.pnt = pnt
        self.shf = shf
        self.mt = mt
        self.pa = pa
        self.pb = pb
    
    @property
    def tab(self):
        return [self.ma, self.mb, self.za, self.zb, self.ia, self.ib, self.q, 
                self.pnt, self.shf, self.mt, self.pa, self.pb]
    


class J_group(list):
    
    def __init__(self, endf):
        self.aj, self.pj, self.kbk, self.kps, nchX6, nch, tab = endf.read_list()
        if nch > 5:
            logging.error("MF2 MT151 : SANDY cannot handle R-matrix limited format with 'nch={}'".format(nch))
            sys.exit()
        for i in range(nch):
            sec = Channel(*tab[i*6:6+i*6])
            self.append(sec)
        _, _, _, self.nrs, nxX6, nx, tab = endf.read_list()
        tab = np.reshape(tab, (nx,6)).T
        er = tab[0]
        gams = tab[1:]
        # Add widths to channels
        for i in range(nch):
            self[i].res = RMLres(er, gams[i])
    
    @property
    def aj(self):
        r"""
        Floating point value of J (spin); sign indicates parity.
        """
        return self.parity * self.j
    
    @aj.setter
    def aj(self, aj):
        self._aj = aj
    
    @property
    def j(self):
        r"""
        Spin.
        """
        return abs(self._aj)
    
    @property
    def parity(self):
        r"""
        Sign of the parity.
        """
        return np.sign(self._aj)

    @property
    def kbk(self):
        """
        Non-zero if background R-matrix exist.
        """
        return self._kbk

    @kbk.setter
    def kbk(self, kbk):
        if kbk != 0:
            logging.error("MF2 MT151 : SANDY cannot handle R-matrix limited format with 'kbk={}'".format(kbk))
            sys.exit()
        self._kbk = kbk

    @property
    def kps(self):
        """
        Non-zero if non-hard-sphere phase shift are to be specified.
        """
        return self._kps

    @kps.setter
    def kps(self, kps):
        if kps != 0:
            logging.error("MF2 MT151 : SANDY cannot handle R-matrix limited format with 'kps={}'".format(kps))
            sys.exit()
        self._kps = kps

    def write(self, ismp=None):
        from sandy.endf import ENDF
        text = []
        tab = np.ravel([ ch.tab for ch in self ])
        text += ENDF.write_list(self.aj, self.pj, self.kbk, self.kps, len(self), tab)
        tab = np.zeros((self.nrs, 6))
        for i,channel in enumerate(self):
            tab[:,0] = channel.res.er
            tab[:,i+1] = channel.res.gam
        tab = tab.reshape(-1)
        text += ENDF.write_list(0, 0, 0, self.nrs, self.nrs, tab)
        return text



class Channel:
    r"""
    RML channel section for a given orbital angular momentum `l`.
    """
    
    def __repr__(self):
        return "<Channel l={}>".format(self.l)
    
    def __init__(self, ppi, l, sch, bnd, ape, apt):
        self.ppi = ppi
        self.l = l
        self.sch = sch
        self.bnd = bnd
        self.ape = ape
        self.apt = apt
    
    @property
    def tab(self):
        r"""
        Channel parameters into a list.
        """
        return [self.ppi, self.l, self.sch, self.bnd, self.ape, self.apt]




class RMLres:
    r"""
    R-matrix limited format resonance.
    """

    def __repr__(self):
        return "<RML resonance>"

    def __init__(self, er, gam):
        self.er = np.array(er)
        self.gam = np.array(gam)
    
#    def tolist(self, ismp=None):
#        if ismp is not None:
#            LIST = self.er[:,ismp].tolist()
#            LIST += self.gam[:,ismp].tolist()
#        else:
#            LIST = self.er.tolist() + self.gam.tolist()
#        return LIST



class MF(MFstandard):
    r"""
    ``MF2`` section of an ``ENDF-6`` file.
    All the corresponding ``MT`` subsections are stored are elements of the 
    dictionary.
    """
    
    def __init__(self, endf):
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Processing ENDF-6 section MF2 : RESONANCE PARAMETERS "))
        logging.debug(80*"+")
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend
    


class MT(list, MTstandard):
    """
    ``MT`` section of file ``MF2``.
    This section is used to store resonance parameters.
    """

    @property
    def nis(self):
        """
        Number of isotopes in the material.
        ``SANDY`` can work only with ``nis=1``.
        """
        return self._nis
    
    @nis.setter
    def nis(self, nis):
        if nis > 1:
            logging.error(self.prefix + "{} isotopes for MAT={}. SANDY can work only with 'nis=1".format(nis, self.mat))
            sys.exit()
        self._nis = nis
    
    @property
    def zai(self):
        r"""
        ``(Z,A)`` designation for an isotope.
        """
        return self._zai

    @zai.setter
    def zai(self, zai):
        if zai != self.za:
            logging.error(self.prefix + "SANDY can work only with one isotope, hence 'zai' must be equal to 'za'")
            sys.exit()
        self._zai = zai

    @property
    def abn(self):
        r"""
        Abundance of an isotope in the material.
        This is a number fraction, not a weight fraction, nor a percent.
        """
        return self._abn

    @abn.setter
    def abn(self, abn):
        if abn != 1:
            logging.error(self.prefix + "SANDY can work only with one isotope, hence 'abn' must be equal to 1.0")
            sys.exit()
        self._abn = abn
    
    @property
    def lfw(self):
        r"""
        Flag indicating whether average fission widths are given in the 
        unresolved resonance region for this isotope.
        - :``lfw=0``: :
            average fission widths are not given
        - :``lfw=1``: :
            average fission widths are given
        """
        return self._lfw
    
    @lfw.setter
    def lfw(self, lfw):
        if lfw == 0:
            logging.debug(self.prefix + "Average fission widths are not given in the unresolved resonance region")
        elif lfw == 1:
            logging.debug(self.prefix + "Average fission widths are given in the unresolved resonance region")
        else:
            logging.error(self.prefix + "'lfw' flag to indicate whether average fission widths are given in the URR must be 0 or 1: found lfw={}".format(lfw))
            sys.exit()
        self._lfw = lfw
    
    @property
    def ner(self):
        """
        Number of resonance energy ranges for isotope.
        """
        return self._ner
        
    
    @ner.setter
    def ner(self, ner):
        if ner < 1:
            logging.error(self.prefix + "At least one resonance energy range must be given: found 'ner={}'".format(ner))
            sys.exit()
        logging.debug(self.prefix + "'ner={}' resonance regions".format(ner))
        self._ner = ner

    @property
    def lru(self):
        r"""
        Flag indicating whether this energy range contains data for resolved 
        or unresolved resonance parameters:
            - :``lru=0``: :
                only the scattering radius is given (``lrf=0``, ``nls=0``, ``lfw=0``)
            - :``lru=1``: :
                resolved resonance parameters are given
            - :``lru=2``: :
                unresolved resonance parameters are given
        """
        return self._lru

    @lru.setter
    def lru(self, lru):
        if lru not in [0,1,2]:
            logging.error(self.prefix + "Format 'lru={}' is not allowed".format(lru))
            sys.exit()
        if lru == 0 and self.ner != 1:
            logging.error(self.prefix + "When only the scattering radius is given ('lru=0'), the number of energy ranges must be 0, not 'ner={}'".format(self.ner))
            sys.exit()
        self._lru = lru

    @property
    def lrf(self):
        r"""
        Flag indicating which representation has been used for the energy 
        range.
        The definition of ``lrf`` depends on the value of ``lru``.
        
        If ``lru=1`` (resolved parameters), then:
            - :``lrf=1``: :
                single-level Breit-Wigner (SLBW)
            - :``lrf=2``: :
                multilevel Breit-Wigner (MLBW)
            - :``lrf=3``: :
                Reich-Moore (RM)
            - :``lrf=4``: :
                Adler-Adler (AA) (Not covered by ``SANDY``)
            - :``lrf=5``: :
                no longer available
            - :``lrf=6``: : 
                no longer available
            - :``lrf=7``: :
                R-Matrix Limited (RML)

        If ``lru=2`` (unresolved parameters), then:
            - :``lrf=1``: :
                only average fission widths are energy-dependent
            - :``lrf=2``: :
                average level spacing, competitive reaction widths, reduced
                neutron widths, radiation widths, and fission widths are 
                energy-dependent.
        """
        return self._lrf

    @lrf.setter
    def lrf(self, lrf):
        if self.lru == 0 and lrf != 0:
            logging.error(self.prefix + "When flag 'lru=0', then 'lrf=0', not 'lrf={}'".format(lrf))
            sys.exit()
        if self.lru == 1 :
            if lrf not in range(8):
                logging.error(self.prefix + "Format 'lrf={}' is not allowed for RRR parameters".format(lrf))
                sys.exit()
            if lrf not in [0,1,2,3,7]:
                logging.error(self.prefix + "SANDY cannot process RRR parameters in format 'lrf={}'".format(lrf))
                sys.exit()
        if self.lru == 2 and lrf not in [1,2]:
            logging.error(self.prefix + "Format 'lrf={}' is not allowed for URR parameters".format(lrf))
            sys.exit()
        self._lrf = lrf

    RRR_formalisms = { 0 : RRR,
                       1 : BW,
                       2 : BW,
                       3 : RM,
                       7 : RML}

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.za, self.awr, self.nis = endf.read_cont([2,3,5])
        self.zai, self.abn, self.lfw, self.ner = endf.read_cont([2,5])
        for i in range(self.ner):
            self.lru = endf.line.CONT.l1
            self.lrf = endf.line.CONT.l2
            if self.lru == 0 or self.lru == 1: # RESOLVED RESONANCES
                sec = self.RRR_formalisms[self.lrf](endf)
                self.append(sec)
            elif self.lru == 2: # UNRESOLVED RESONANCES
                nurr = self.ner - len(self) # Number of URR sections
                logging.debug("MF{} MT{} URR : Found {} unresolved resonance regions".format(self.mf, self.mt, nurr))
                if self.mf == 32:
                    logging.debug("MF{} MT{} URR : SANDY cannot process covariance section for the unresolved resonance region".format(self.mf, self.mt))
                self.urr = URR.read_endf(endf)

    def write(self, ismp=None):
        r"""
        Return the text of the ``MT`` section as a list of strings.
        
        Inputs:
            - ismp :
                (integer) number of the sample

        Ouputs:
            - :``text``: :
                list of strings
        """
        from sandy.endf import ENDF
        text = []
        text += ENDF.write_cont(self.za, self.awr, 0, 0, self.nis, 0)
        text += ENDF.write_cont(self.zai, self.abn, 0, self.lfw, self.ner, 0)
        for section in self:
            text += section.write(ismp=ismp)
        if hasattr(self, 'urr'):
            text += self.urr
        text = ENDF.add_control(text, self.mat, self.mf, self.mt)
        return text
    
    def sample_resonances(self, mt32):
        r"""
        Loop the covariance sections ``mt32`` and draw samples of the 
        resonance parameters.
        Then, for each resonance with samples, look for its analogous in the 
        ``self`` subsections and replace it.
        
        Print a warning if the covariance is given for a resonance that is 
        not present in ``self`` (MF2).
        
        Inputs:
            - :``mt32``: :
                (mf2.MT instance) covariance sections
        """
        from sandy.sandy_input import options
        nsmp = options['samples']
        for sec in mt32:
            if hasattr(sec, "covsec"):
                # Draw samples from the covariance section
                sec.covsec.sampling(nsmp)
                sec.covsec.smp[0].tolist()
                for res in sec.covsec.smp:
                    found = [ erange.replace_resonance(res) for erange in self ]
                    if not np.any(found):
                        logging.debug(self.prefix + "Covariance given for resonance (E={:.5e}, J={}) but resonance was not found in MF2".format(res.er, res.aj))
    #                        if "stdmax" in options:
    #                            tolerance = options["stdmax"]
    #                            erange.ap.filter_samples(tolerance, sec.covsr)

    def sample_sr(self, mt32):
        r"""
        Loop the covariance sections ``mt32`` and draw samples of the 
        scattering radius.
        The samples are applied to all the scattering radius values given in 
        the ``self`` subsection, independently from the energy range or 
        *l*-value.
        
        Inputs:
            - :``mt32``: :
                (``mf2.MT`` instance) covariance sections
        """
        from sandy.sandy_input import options
        nsmp = options['samples']
        for sec in mt32:
            if hasattr(sec, "covsr"):
                smp = np.ravel(sec.covsr.sampling(nsmp))
                for erange in self:
                    if isinstance(erange, BW) or isinstance(erange, RM):
                        logging.debug(erange.prefix + "Absolute uncertainty 'dap={}' given to scattering radius 'ap={}'".format(sec.covsr.std[0], erange.ap))
                    erange.ap.smp = smp
                    erange.ap.replace_negative()
                    if isinstance(erange, RM):
                        for lsec in erange:
                            logging.debug(lsec.prefix + "Absolute uncertainty 'dap={}' given to scattering radius 'apl={}'".format(sec.covsr.std[0], lsec.apl))
                            lsec.apl.smp = smp
                            lsec.apl.replace_negative()
                return
                
