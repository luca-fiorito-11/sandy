# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:31:01 2017

@author: lfiorito
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 11:11:14 2015

@author: lfiorito
"""
import numpy as np
from sandy.records import Tab1, Tab2, Samples
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard
import logging


class LPC(Tab1):
    
    mu = np.linspace(-1, 1, 101)
    
    @classmethod
    def read_endf(cls, endf):
        r"""
        Initialize an instance given only the grid size.
        Define an unitary grid that starts from zero as :math:`x` and a zero 
        vector as :math:`f(x)`.
        
        .. Tip::
            x = 0  1  2  3  ... N
            y = 0. 0. 0. 0. ... 0.

        Use zero interpolation over the whole grid.
        
        Inputs:
            - N :
                (integer) size of the grid
        """
        nm = 64
        t, e, lt, _, nl, _, a_coeff = endf.read_list()
        xx = np.arange(1, nm+1)
        yy = np.zeros(nm)
        yy[:nl] = a_coeff
        inst = cls(xx, yy, nbt=nm, intscheme=1)
        inst.nl = nl
        inst.t = t
        inst.lt = lt
        return inst
    
    @property
    def nl(self):
        r"""
        Highest order of Legendre polynomials for this representation.
        If ``nl`` is not set, use the maximum value :math:`nl=64`.
        """
        if hasattr(self, "_nl"):
            return self._nl
        else:
            return 64

    @nl.setter
    def nl(self, nl):
        self._nl = nl
    
    @property
    def t(self):
        r"""
        Temperature (in K), this value is normally zero.
        If ``t`` is not set, use :math:`t=0`.
        """
        if hasattr(self, "_t"):
            return self._t
        else:
            return 0

    @t.setter
    def t(self, t):
        self._t = t

    @property
    def lt(self):
        r"""
        Test for temperature dependence, this value is normally zero.
        If ``lt`` is not set, use :math:`lt=0`.
        """
        if hasattr(self, "_lt"):
            return self._lt
        else:
            return 0

    @lt.setter
    def lt(self, lt):
        self._lt = lt        

    @property
    def coefficients(self):
        r"""
        Get coefficients :math:`w_l` (with :math:`l \in [0, NL]`) of the 
        Legendre polynomial series representing the angular distribution

        .. math::
            
            f(\mu,E) = \sum_{l=0}^{NL} w_l(E) P_l(\mu)
        
        with
        
        .. math::
            
            w_l(E) = \frac{2l+1}{2} a_l(E)
        """
        shape = (1,self.shape[1]) if self.ndim == 2 else (1,)
        aa = np.concatenate((np.ones(shape), self.yy[:self.nl]))  # add a0 = 1.0
        return (aa.T * ((2*np.arange(self.nl+1) + 1)/2.)).T
    
    @property
    def tpd(self):
        r"""
        Tabulated probability distribution :math:`f(\mu,E_i)` for a given 
        energy :math:`E_i`.
        The distribution is obtained from the Legendre polynomial coefficients 
        as
        
        .. math::

            f(\mu,E_i) = \sum_{l=0}^{NL} \frac{2l+1}{2} a_l(E_i) P_l(\mu)
        """
        return np.polynomial.legendre.legval(self.mu, self.coefficients).T

    def write_endf(self, c1=0, c2=0, l1=0, l2=0, ismp=None):
        r"""
        Write ``Tab1`` record for angular distribution expressed as Legendre 
        polynomials according to ``ENDF-6`` format.
        """
        from sandy.endf import ENDF
        yy = self if ismp is None else self[:,ismp]
        if self.xx.size != yy.size:
            raise NotImplementedError("Legendre Polinomial coefficients have ",
                                      "the wrong size")
        yy = yy[:self.nl]
        text = ENDF.write_list(0, c2, 0, 0, 0, yy)
        return text

    def plot(self, ax, color='b', label=None, **kwargs):
        r"""
        Add plot of the :math:`f(\mu,E_i)` distribution.
        """
        from sandy.plot import plot1d
        plot1d(ax, self.mu, self.tpd, color=color, label=label)

    def fig(self):
        r"""
        Produce ``matplotlib`` figure of the :math:`f(\mu,E_i)` distribution.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        self.plot(ax)
        ax.set_xlim([-1, 1])
        ax.set_xscale("linear")
        ax.set_yscale("linear")
        ax.grid()
        ax.set_xlabel(r"$\mu$")
        ax.set_ylabel(r"$f(\mu,E)$")
        fig.show()



class LPCSmp(LPC, Samples):
    pass



class Tab2_LPC(Tab2):
    """
    Subclass of class ``records.Tab2``.
    
    Methods ``records.Tab2.read_endf`` and ``records.Tab2.write_endf`` 
    are redefined because the tabulated function :math:`f(x,y)` along the 
    second dimension, that is :math:`f(x_1,y), f(x_2,y), \dots`, 
    is provided as a ``List`` record and not a ``Tab1`` record.
    """
    
    @property
    def prefix(self):
        return "MF{} MT{} LPC : ".format(self.mf, self.mt)
    

    @property
    def nm(self):
        r"""
        Maximum order of Legendre polynomial that is required to describe 
        the angular distributions in either the center-of-mass or the 
        laboratory system.
        ``nm`` should be an even number.
        """
        if hasattr(self, "_nm"):
            return self._nm
        else:
            return 64
    
    @nm.setter
    def nm(self, nm):
        if nm == 0:
            nm = 64
        self._nm = nm

    def legendre_coefficient(self, l):
        r"""
        Extract a single energy-dependent Legendre polynomial coefficient.
        
        Inputs:
            - :``l``: :
                (scalar integer) index of the Legendre polynomial coefficient
        
        Outputs:
            - :``inst``: :
                (`mf4.P` instance) energy-dependent Legendre polynomial coefficient
        """
        try:
            l = int(abs(l))
        except (ValueError, TypeError):
            raise NotImplementedError("Legendre polynomal index must be an integer")
        if l > self.nm:
            raise NotImplementedError("Requested order of the Legendre ",
                                      "polynomial exceeds the maximum allowed")
        yy = [ y[l-1] for y in self ]
        mytype = Samples if np.ndim(yy) == 2 else Tab1 
        inst = mytype(self.xx, yy, nbt=self.nbt, intscheme=self.intscheme)
        return inst

    def perturb_lpc(self, pert, l):
        r"""
        Given energy-dependent relative samples for a Legendre Polynomial 
        coefficient perturb the best estimate value.
        
        .. Important::
            the samples are assumed to be in relative unit.
        
        Inputs:
            - :``pert``: :
                (Samples instance) energy-dependent relative samples for the 
                Legendre Polynomial coefficient
            - :``l``: :
                (scalar integer) index of the Legendre polynomial coefficient
                (:math:`l > 0`)
        """
        for i,p in enumerate(pert):
            self[i][l-1] *= (1 + pert[i])

    def fig(self, xx):
        r"""
        Produce ``matplotlib`` figure of :math:`f(\mu,E)` distribution for 
        given incoming-neutron energy.
        
        Inputs:
            - :``xx``: :
                (scalar or 1d-array) incoming-neutron's single energy or array 
                of energies
        """
        import matplotlib.pyplot as plt
        from sandy.plot import colors
        xx = np.atleast_1d(xx).astype(float)
        yy = self(xx)
        fig, ax = plt.subplots()
        for i,(x,y) in enumerate(zip(xx,yy)):
            y.plot(ax, color=colors[i], label=r"$E={:.5e} eV$".format(x))
        ax.grid()
        ax.set_xlabel(r"$\mu$")
        ax.set_ylabel(r"$f(\mu,E)$")
        ax.set_xlim([-1, 1])
        ax.legend(loc='best', fancybox=True, shadow=True)
        fig.show()



class MF(MFstandard):
    """
    ``MF4`` section of an ``ENDF-6`` file.
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
    
    def sections_with_attr(self, attr):
        r"""
        Set of ``MT`` values of the sections with a requested attribute.
        """
        return set([ mt for mt in self if hasattr(self[mt], attr) ])

    def add_samples(self, mf34):
        r"""
        """
        from sandy.sandy_input import options
        nsmp = options['samples']
        extra_points = mf34.union_grid.add_extra_points()
        for mt in sorted(self.sections_with_attr('lpc') & mf34.keys()):
            self[mt].lpc = self[mt].lpc.add_points(extra_points)
            lpcsmp = self[mt].lpc.init_smp(nsmp, LPCSmp)
            self[mt].cov = {}
            for ((x, l), (xi, li)), sec34 in mf34.sorted_items():
                if hasattr(sec34, "smp"):
                    rel_pert = sec34.smp.reinit(lpcsmp.xx)
                    lpcsmp.perturb_lpc(rel_pert, l)
                    self[mt].cov[l] = sec34
            # wait that all Ls have been perturbed
            self[mt].smp = lpcsmp
            
    def plot_samples(self, filename):
        r"""
        Plot samples standard deviation against ``MF34`` standard deviation 
        of the first y-axis.
        Plot the cross section on the second y-axis.
        """
        import matplotlib.pyplot as plt
        from sandy.sandy_input import outdir, options
        from sandy.plot import save
        from os.path import join
        nsmp = options['samples']
        directory = join(outdir, "plots-" + filename)
        for mt in self.sections_with_attr('smp'):
            for l in self[mt].cov:
                fig, ax = plt.subplots()
                std = self[mt].cov[l].Std
                std *= 100
                std.plot(ax, step=True, color='r', label='ENDF-6')
                LegCoeff = self[mt].smp.legendre_coefficient(l=l)
                rstd = LegCoeff.Rstd
                rstd *= 100
                rstd.plot(ax, step=True, color='b', label='samples')
                ax.set_title(r"MF4 MT{}, P{}, {} samples".format(mt, l, nsmp))
                ax.set_xscale('log')
                ax.axes.set_ylabel('stdev (%)')
                ax.axes.set_xlabel('energy (eV)')
                ymin, ymax = ax.get_ylim()
                ax.set_ylim([0, min(ymax,100)])
                ax.grid()
                ax.legend(loc='best', fancybox=True, shadow=True)
                ax2 = ax.twinx()
                LegCoeff.Mean.plot(ax2, alpha=0.5, color='0.75')
                ax2.set_ylabel('P{} mean (-)'.format(l), color='0.75')
                ax2.tick_params('y', colors='0.75')
                ax2.set_xscale('log')
                ax2.set_yscale('linear')
                figname = join(directory, "samples_mf4_mt{}_p{}".format(mt, l))
                save(figname, fig, 'pdf')



class MT(MTstandard):
    r"""
    ``MT`` section of file ``MF34``.
    This section contains the angular distributions of secondary 
    particles.
    
    The structure of a section depends on the ``ltt`` value, which defines 
    the representation used for the angular distributions.
    """

    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.za, self.awr, self.ltt = endf.read_cont([2,4,5])
        self.awr, self.li, self.lct, nm = endf.read_cont([0,4])
        if self.ltt == 0: # Isotropic distribution
            # Used for ENDF/B-VII.1 U-235 MT18
            pass
        if self.ltt == 1 or self.ltt == 3: # Legendre polynomial
            self.lpc = Tab2_LPC.read_endf(endf, LPC)
            self.lpc.nm = nm
        if self.ltt == 2 or self.ltt == 3: # Tabulated distribution
            self.tpd = Tab2.read_endf(endf)
    
    @property
    def ltt(self):
        r"""
        Flag to specify the representation used and it may have the 
        following values.
            - :``ltt=0``: :
                all angular distributions are isotropic
            - :``ltt=1``: :
                the data are given as Legendre expansion coefficients, 
                :math:`a_l (E)`
            - :``ltt=2``: :
                the data are given as tabulated probability distributions, 
                :math:`f(\mu, E)`
            - :``ltt=3``: :
                low energy region is represented by as Legendre coefficients; 
                higher region is represented by tabulated data.
        """
        return self._ltt
    
    @ltt.setter
    def ltt(self, ltt):
        self._ltt = ltt
    
    @property
    def lct(self):
        r"""
        Flag to specify the frame of reference used.
            - :``lct=1``: :
                the data are given in the laboratory system
            - :``lct=2``: :
                the data are given in the center-of-mass system
        """
        return self._lct
    
    @lct.setter
    def lct(self, lct):
        self._lct = lct
    
    @property
    def li(self):
        """
        Flag to specify whether all the angular distributions are isotropic.
            - :``li=0``: :
                not all isotropic
            - :``li=1``: :
                all isotropic
        """
        return self._li
    
    @li.setter
    def li(self, li):
        self._li = li
    
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
        text += ENDF.write_cont(self.za, self.awr, 0, self.ltt, 0, 0)
        nm = self.lpc.nm if self.ltt == 3 else 0
        text += ENDF.write_cont(0, self.awr, self.li, self.lct, 0, nm)
        if self.ltt == 1 or self.ltt == 3:
            if ismp is not None and hasattr(self, 'smp'):
                text += self.smp.write_endf(ismp=ismp)
            else:
                text += self.lpc.write_endf()
        if self.ltt == 2 or self.ltt == 3:
            text += self.tpd.write_endf()
        text = ENDF.add_control(text, self.mat, self.mf, self.mt)
        return text