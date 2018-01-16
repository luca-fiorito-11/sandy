# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:16:23 2017

@author: lfiorito
"""
from sandy.records import Tab1, Tab2, Samples
import numpy as np
import logging
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard


class TED(Tab1):
    r"""
    Tabulated energy distribution.
    """
    
    @property
    def integral(self):
        r"""
        Integral :math:`I` of the tabulated energy distribution
        .. math::
            I = \int_{E_{out}^{min}}^{E_{out}^{max}} g(E\rightarrow E_{out}) dE_{out} 
        """
        fx = (self[1:] + self[:-1])/2.
        return fx.T.dot(self.dx)

    def normalize(self, norm=1.):
        r"""
        Normalize Tabulated energy distribution.

        ..Important::
            this function works for 1d- and 2d-arrays.
        
        Inputs:
            - norm :
                (scalar float) normalization constant, default is 1.0
        """
        self *= norm/self.integral
    
    def is_normalized(self, norm=1., rtol=1e-5):
        return np.isclose(self.integral(), norm, rtol=rtol)
    


class TEDSmp(TED, Samples):
    pass



class ATF(Tab2):
    r"""
    Arbitrary Tabulated Function (LF=1).
    
     - 1st dimension :
         incoming neutron energy :math:`E_{in}`
     - 2nd dimension :
         outgoing neutron energy :math:`E_{out}`
     - function :
         partial energy distribution :math:`g(E_{in}\rightarrow E_{out})
     
    The partial energy distribution is integrated over each energy bin and must 
    be normalized to one.
    According to the ``ENDF-6`` manual, the interpolation scheme used between 
    incident energy points, and between secondary energy points, should be 
    linear-linear.
    
    .. Important::
        Note that the incident energy mesh for :math:`p_k(E)` does not have to 
        be the same as the mesh used to specify the energy distributions.
    """
    
    @classmethod
    def read_endf(cls, endf):
        return super(ATF, cls).read_endf(endf, TYPE=TED)

    @property
    def lf(self):
        r"""
        Flag specifying the energy distribution law used for a particular 
        subsection (partial energy distribution).
        """
        return 1

    @property
    def p(self):
        r"""
        :math:`p(E_{in})` is the fractional probability that this distribution can 
        be used at energy :math:`E_{in}`.
        """
        return self._p

    @p.setter
    def p(self, p):
        self._p = p

    def write_endf(self, **kwargs):
        r"""
        Write function in ``ENDF-6`` format
        
        Outputs:
            - :``text``: :
                (list of strings)
        """
        text = []
        text += self.p.write_endf(l2=self.lf)
        text += super().write_endf(**kwargs)
        return text



class GES:
    r"""
    General Evaporation Spectrum ``(LF=5)``.
    
    Found in:
        - ``ENDF/B-VII.1`` Cf252
        - ``ENDF/B-VII.1`` Th227
    """

    @classmethod
    def read_endf(cls, endf):
        theta = Tab1.read_endf(endf)
        g = Tab1.read_endf(endf)
        return cls(theta, g)
    
    def __init__(self, theta, g):
        self.theta = theta
        self.g = g

    @property
    def lf(self):
        r"""
        Flag specifying the energy distribution law used for a particular 
        subsection (partial energy distribution).
        """
        return 5

    @property
    def u(self):
        return self._u
    
    @u.setter
    def u(self, u):
        self._u = u

    @property
    def p(self):
        r"""
        :math:`p(E_{in})` is the fractional probability that this distribution 
        can be used at energy :math:`E_{in}`.
        """
        return self._p

    @p.setter
    def p(self, p):
        self._p = p

    def write_endf(self, **kwargs):
        r"""
        Write function in ``ENDF-6`` format
        
        Outputs:
            - :``text``: :
                (list of strings)
        """
        text = []
        text += self.p.write_endf(c1=self.u, l2=self.lf)
        text += self.theta.write_endf()
        text += self.g.write_endf()
        return text



class SMFS:
    """
    Simple Maxwellian Fission Spectrum ``(LF=7)``.

    Found in:
        - ``ENDF/B-VII.1`` Ra-223
    """
    
    @classmethod
    def read_endf(cls, endf):
        theta = Tab1.read_endf(endf)
        return cls(theta)
    
    def __init__(self, theta):
        self.theta = theta

    @property
    def lf(self):
        r"""
        Flag specifying the energy distribution law used for a particular 
        subsection (partial energy distribution).
        """
        return 7

    @property
    def u(self):
        return self._u
    
    @u.setter
    def u(self, u):
        self._u = u

    @property
    def p(self):
        r"""
        :math:`p(E_{in})` is the fractional probability that this distribution can 
        be used at energy :math:`E_{in}`.
        """
        return self._p

    @p.setter
    def p(self, p):
        self._p = p

    def write_endf(self, **kwargs):
        r"""
        Write function in ``ENDF-6`` format
        
        Outputs:
            - :``text``: :
                (list of strings)
        """
        text = []
        text += self.p.write_endf(c1=self.u, l2=self.lf)
        text += self.theta.write_endf()
        return text



class ES:
    """
    Evaporation Spectrum ``(LF=9)``.
    
    Found in:
        - ``ENDF/B-VII.1`` C-nat
    """

    @classmethod    
    def read_endf(cls, endf):
        theta = Tab1.read_endf(endf)
        return cls(theta)

    def __init__(self, theta):
        self.theta = theta

    @property
    def lf(self):
        r"""
        Flag specifying the energy distribution law used for a particular 
        subsection (partial energy distribution).
        """
        return 9

    @property
    def u(self):
        return self._u
    
    @u.setter
    def u(self, u):
        self._u = u

    @property
    def p(self):
        r"""
        :math:`p(E_{in})` is the fractional probability that this distribution can 
        be used at energy :math:`E_{in}`.
        """
        return self._p

    @p.setter
    def p(self, p):
        self._p = p
    
    def write_endf(self, **kwargs):
        r"""
        Write function in ``ENDF-6`` format
        
        Outputs:
            - :``text``: :
                (list of strings)
        """
        text = []
        text += self.p.write_endf(c1=self.u, l2=self.lf)
        text += self.theta.write_endf()
        return text



class EDWS:
    """
    Energy-Dependent Watt Spectrum ``(LF=11)``.
    
    Found in:
        - ``ENDF/B-VII.1`` U-233
    """
    
    @classmethod    
    def read_endf(cls, endf):
        a = Tab1.read_endf(endf)
        b = Tab1.read_endf(endf)
        return cls(a, b)

    def __init__(self, a, b):
        self.a = a
        self.b = b

    @property
    def lf(self):
        r"""
        Flag specifying the energy distribution law used for a particular 
        subsection (partial energy distribution).
        """
        return 11

    @property
    def u(self):
        return self._u
    
    @u.setter
    def u(self, u):
        self._u = u

    @property
    def p(self):
        r"""
        :math:`p(E_{in})` is the fractional probability that this distribution can 
        be used at energy :math:`E_{in}`.
        """
        return self._p

    @p.setter
    def p(self, p):
        self._p = p

    def write_endf(self, **kwargs):
        r"""
        Write function in ``ENDF-6`` format
        
        Outputs:
            - :``text``: :
                (list of strings)
        """
        text = []
        text += self.p.write_endf(c1=self.u, l2=self.lf)
        text += self.a.write_endf()
        text += self.b.write_endf()
        return text


class MNFS:
    """
    Energy-Dependent Fission Neutron Spectrum (Madland and Nix) ``(LF=12)``.
    """

    @classmethod    
    def read_endf(cls, endf):
        efl = endf.line.CONT.c1
        efh = endf.line.CONT.c2
        Tm = Tab1.read_endf(endf)
        return cls(efl, efh, Tm)
    
    def __init__(self, efl, efh, Tm):
        self.efl = efl
        self.efh = efh
        self.Tm = Tm

    @property
    def lf(self):
        r"""
        Flag specifying the energy distribution law used for a particular 
        subsection (partial energy distribution).
        """
        return 12

    @property
    def p(self):
        r"""
        :math:`p(E_{in})` is the fractional probability that this distribution can 
        be used at energy :math:`E_{in}`.
        """
        return self._p

    @p.setter
    def p(self, p):
        self._p = p

    def write_endf(self, **kwargs):
        r"""
        Write function in ``ENDF-6`` format
        
        Outputs:
            - :``text``: :
                (list of strings)
        """
        text = []
        text += self.p.write_endf(l2=self.lf)
        text += self.Tm.write_endf(c1=self.efl, c2=self.efh)
        return text



class MF(MFstandard):
    """
    ``MF5`` section of an ``ENDF-6`` file.
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
    
    def add_samples(self, mf35):
        from sandy.sandy_input import options
        nsmp = options['samples']
        for mt in sorted(mf35.keys() & self.keys()):
            for sec in self[mt].list_ATF: # List only the arbitrary tabulated functions
                sec.chi = sec.chi.add_points(mf35[mt].xx, sec.prefix)
                atfsmp = sec.chi.init_smp(nsmp, TEDSmp)
                for i,(e_in,chi) in enumerate(atfsmp.items()):
                    try:
                        # Do not change SMP, it can be reused for other energies
                        smp = mf35[mt].get_samples_by_energy(e_in)
                        cov = mf35[mt].get_cov_by_energy(e_in)
                    except AttributeError:
                        logging.debug("Samples were not found for incoming energy E={:.5e} eV".format(e_in))
                        continue
                    # Change samples and energy distribution secondary grid
                    ux = smp.xx.add_extra_points(sec.chi[i].xx)
                    sec.chi[i] = sec.chi[i].reinit(ux)
                    chi = chi.reinit(ux)
                    chi += smp.reinit(ux)
                    chi.replace_negative()
                    chi.normalize()
                    # Add covariance matrix, this is not really a preatty way
                    chi.cov = cov
                    atfsmp[i] = chi
                sec.smp = atfsmp
                if 'stdmax' in options:
                    sec.filter_samples(options['stdmax'])

    def plot_samples(self, filename, **options):
        """
        Call the samples-plotting method of each `MT` section.
        
        Inputs:
            - filename :
                (string) name of the `ENDF-6` file
        """
        for mt in self:
            for section in self[mt]:
                section.fig_samples(filename)



class MT(list, MTstandard):
    """
    ``MT`` section of ``MF5``.
    This section is used to store subsections (list) with partial 
    energy-distributions of secondary particles.
    """
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.za, self.awr, nk = endf.read_cont([2,3,5])
        for k in range(nk):
            self.append(SubSection(endf))
    
    @property
    def list_ATF(self):
        r"""
        List of the ``mf5.SubSections`` with attribute ``chi`` of type 
        ``mf5.ATF``.
        """
        return [ item for item in self if type(item.chi) is ATF ]
                
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
        text += ENDF.write_cont(self.za, self.awr, 0, 0, len(self), 0)
        for section in self:
            if ismp is not None and hasattr(section, "smp"):
                text += section.smp.write_endf(ismp=ismp)
            else:
                text += section.chi.write_endf()
        text = ENDF.add_control(text, self.mat, self.mf, self.mt)
        return text



class SubSection:
    r"""
    SUbsection of ``MT`` section of ``MF5``.
    This section is used to store individual partial energy-distributions of 
    secondary particles.
    """
    laws = { 1 : ATF,
             5 : GES,
             7 : SMFS,
             9 : ES,
            11 : EDWS,
            12 : MNFS}
    
    def __repr__(self):
        return "<MF{} MT{} CHI>".format(self.mf, self.mt)
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        u = endf.line.CONT.c1
        lf = endf.line.CONT.l2
        p = Tab1.read_endf(endf)
        chi = self.laws[lf].read_endf(endf)
        chi.p = p
        if lf in [5, 7, 9, 11, 12]:
            chi.u = u
        self.chi = chi
    
    @property
    def prefix(self):
        return "MF{} MT{} CHI : ".format(self.mf, self.mt)

    def filter_samples(self, threshold):
        r"""
        All samples that have the standard deviation above a given threshold 
        are canceled and must be replaced with the mean.
        
        Inputs:
            - :``threshold``: :
                (scalar) threshold of relative standard deviation
        """
        if hasattr(self, 'smp'):
            for e_in,chi in self.smp.items():
                emask, smask = chi.stdmax(threshold, rstd=True)
                mask = np.in1d(chi.xx, emask)
                i = self.chi.where(e_in)
                mean = self.chi[i].interp(chi.xx)
                chi[mask] = mean[mask, np.newaxis]
                if len(emask) > 0:
                    logging.info("MF5 MT{} E_in={:.3e} eV : relative stdev > {} at the following energy points:".format(self.mt, e_in, threshold))
                for x,s in zip(emask, smask):
                    logging.debug(" - E = {:5e} eV --> rstd = {:.2f} %".format(x, s*100.))
    
    def fig_samples(self, filename=None):
        r"""
        Plot samples standard deviation against ``MF35`` standard deviation 
        of the first y-axis.
        Plot the energy-distribution also on the first y-axis.
        """
        from sandy.functions import div0
        import matplotlib.pyplot as plt
        from os.path import join
        from sandy.plot import save
        from sandy.sandy_input import outdir
        if not hasattr(self, 'smp'):
            return
        for i,(e_in,smp) in enumerate(self.smp.items()):
            if not hasattr(smp, "cov"):
                continue
            xmin = self.chi[i].xx[self.chi[i].xx>0][0]
            xmax = self.chi[i].xx[-1]
            xx = np.logspace(np.log10(xmin),np.log10(xmax),1000)
            mean = self.chi[i].reinit(xx)
            fig, ax = plt.subplots()
            # Covariance
            std = smp.cov.Std.reinit(xx)
            rstd = Tab1(xx, div0(std, mean)*100.)
            rstd.plot(ax, step=True, label='ENDF-6', color='r')
            # Samples
            stds = smp.Std.reinit(xx)
            rstds = Tab1(xx, div0(stds, mean)*100.)
            rstds.plot(ax, step=True, label='samples')
            ax.set_title(r"MF5 MT{}, E={:.3e} eV, {} samples".format(self.mt, e_in, smp.nsmp))
            ax.set_ylabel('stdev (-)')
            ax.set_xlabel('energy (eV)')
            ax.set_xlim([xmin, xmax])
            ymin, ymax = ax.get_ylim()
            ax.set_ylim([0, min(ymax,100)])
            ax.grid()
            ax.legend(loc='best', fancybox=True, shadow=True)
            ax2 = ax.twinx()
            i = self.chi.where(e_in)
            mean.plot(ax2, alpha=0.5, color='0.75')
            ax2.set_ylabel(r'pdf (-)', color='0.75')
            ax2.tick_params('y', colors='0.75')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            if filename:
                directory = join(outdir, "plots-" + filename)
                figname = join(directory, "samples_mf{}_mt{}_{:.3e}".format(self.mf, self.mt, e_in))
                save(figname, fig, 'pdf')
            else:
                fig.show()