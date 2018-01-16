# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 17:43:20 2016

@author: lfiorito
"""
import logging
import numpy as np
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard


redundant_xs = {1, 3, 4, 18, 27, 101, 103, 104, 105, 106, 107, 452}


class MF(MFstandard):
    r"""
    ``MF3`` section of an ``ENDF-6`` file.
    All the corresponding ``MT`` subsections are stored are elements of the 
    dictionary.
    """
    
    def __init__(self, endf):
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Processing ENDF-6 section MF3 : CROSS SECTIONS "))
        logging.debug(80*"+")
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            if mt == 1:
                self[mt] = MT1(endf)
            elif mt == 3:
                self[mt] = MT3(endf)
            elif mt == 4:
                self[mt] = MT4(endf)
            elif mt == 18:
                self[mt] = MT18(endf)
            elif mt == 27:
                self[mt] = MT27(endf)
            elif mt == 101:
                self[mt] = MT101(endf)
            elif mt == 103:
                self[mt] = MT103(endf)
            elif mt == 104:
                self[mt] = MT104(endf)
            elif mt == 105:
                self[mt] = MT105(endf)
            elif mt == 106:
                self[mt] = MT106(endf)
            elif mt == 107:
                self[mt] = MT107(endf)
            else:
                self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend
        self._set_daughters()

    @property
    def union_grid(self):
        r"""
        Union energy grid for the cross sections in ``MF3``.
        """
        from sandy.functions import union_grid
        return union_grid(*[self[mt].xs.xx for mt in self])

    @property
    def parents(self):
        """
        List of ``MT`` redundant cross sections.
        
        To be a *parent* it is not only necessary to be a redundant cross 
        section, but its components must be present in the ``ENDF-6`` file.
        For example, ``MT4`` (total inelastic cross section) is a *parent* 
        if any ``MT`` bewteen 51 to 91 is present (the inelastic cross sections 
        with the residual in different inelastic states).
        If not, ``MT4`` is not a *parent*.
        """
        parents = []
        for mt,mtsec in self.items():
            if hasattr(mtsec, "components"):
                if len(mtsec.daughters) > 0:
                    parents.append(mt)
        return parents
    
    def _set_daughters(self):
        """
        Assign *daughters* to the redundant cross section.
        For further info, see docstring of attribute 
        ``mf3.MTRedundant.daughters``.
        
        .. Important:
            This method assigns attribute ``daughters`` to the ``mf3.MT`` 
            objects
        """
        for mt in self:
            if hasattr(self[mt], "components"):
#                import pytest
#                pytest.set_trace()
                self[mt].daughters = self.keys()

    def is_unionized(self, err=True):
        """
        Check if all the cross sections are *unionized* (defined on a common 
        union grid).
        
        Inputs:
            - :``err``: :
                (boolean) flag that raises error or not
        
        Outputs:
            - :``is_unionized``: : 
                (boolean) ``True`` if the cross sections are *unionized*, 
                otherwise ``False``
        """
        is_unionized = []
        for mt,mtsec in sorted(self.items()):
            check = (np.in1d(mtsec.xs.xx, self.union_grid)).all()
            if not check:
                raise NotImplementedError("cross section values cannot be summed because 'MT{}' is not unionized".format(mt))
            is_unionized.append(check)
        is_unionized = np.all(is_unionized)
        return is_unionized

    def sum_rules(self):
        r"""
        Reconstruct redundant cross sections as the sum of their components.
        The function works only if all the tabulated functions are unionized.
        """
        self.is_unionized()
        logging.info(self.prefix + "Reconstruct redundant cross sections")
        for mt in sorted(self.parents, reverse=True):
            items = [ self[dau] for dau in self[mt].daughters ]
            self[mt].SUM(*items)
        
    def add_samples(self, mf33):
        r"""
        Add samples from ``MF33`` to ``MF3``.

        For a given ``MT`` section in ``MF33``:
            * assign the samples to the corresponding ``MT`` section in ``MF3``
            * assign the covariance to the corresponding ``MT`` section in ``MF3``
            * if the ``MT`` section refer to a redundant cross section, then ...
            * ... if none of its component has already got samples ...
            * ... assign the samples to the redundant cross section components 
            * ... assign the covariance to the redundant cross section components 
            * eventually, apply summation rules

        .. Important::
            All the samples in ``MF3`` are tabulated on a unionized grid.

        Inputs:
            - :``mf33``: :
                (``mf33.MF`` instance) dictionary containing ``mf33.MT`` 
                instances
        """
        for mt, mtsec in sorted(mf33.items(), reverse=True):
            if mt not in self:
                logging.warn(mtsec.prefix + "Section MF{} MT{} was not found".format(self.mf, mt))
                continue
            smp = mtsec.smp
            cov = mtsec[mt]
            self[mt]._add_samples(smp)
            self[mt].cov = cov
            if hasattr(self[mt], "components"):
                existing_components = self.keys() & self[mt].components
                if not np.any([ hasattr(self[x], "smp") for x in existing_components ]):
                    for x in sorted(existing_components):
                        logging.debug(mtsec.prefix + "Assign samples to MF{} MT{}".format(self.mf, x))
                        self[x]._add_samples(smp)
                        self[x].cov = cov
#            if hasattr(self[mt], "daughters"):
#                # Check that no daughter has samples
#                if not np.any([ hasattr(self[x], "smp") for x in self[mt].daughters ]):
#                    for dau in self[mt].daughters:
#                        logging.debug(mtsec.prefix + "Assign samples to MF{} MT{}".format(self.mf, dau))
#                        self[dau]._add_samples(smp)
#                        self[dau].cov = cov
        self.sum_rules()
        
    def plot_samples(self, filename):
        """
        Call the samples-plotting method of each ``MT`` section.
        
        Inputs:
            - :``filename``: :
                (string) name of the ``ENDF-6`` file
        """
        for mt,mtsec in sorted(self.items()):
            mtsec.fig_samples(filename)
    


class MT(MTstandard):
    """
    ``MT`` section of ``MF3``.
    This section is used to store tabulated neutron-induced cross sections.
    """
    @property
    def qm(self):
        """
        Mass-difference Q value (eV).
        """
        return self._qm
    
    @qm.setter
    def qm(self, qm):
        self._qm = qm

    @property
    def qi(self):
        """
        Reaction Q value for the (lowest energy) state defined by the 
        given `MT` value in a simple two-body reaction or a breakup 
        reaction.
        """
        return self._qi
    
    @qi.setter
    def qi(self, qi):
        self._qi = qi
    
    @property
    def lr(self):
        r"""
        Complex or *breakup* reaction flag, which indicates that additional 
        particles not specified by the ``MT`` number will be emitted.
        """
        return self._lr
    
    @lr.setter
    def lr(self, lr):
        self._lr = lr
    
    def __init__(self, endf):
        from sandy.records import Tab1
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.l2 = endf.read_cont([0,1,2,4,5])
        self.qm = endf.line.CONT.c1
        self.qi = endf.line.CONT.c2
        self.lr = endf.line.CONT.l2
        self.xs = Tab1.read_endf(endf)
        logging.debug(self.prefix + "emin={:.5e} eV, nbt={}, intscheme={}".format(self.xs.xx[0], self.xs.nbt, self.xs.intscheme))
    
    def _add_samples(self, smp):
        r"""
        Perturb cross sections according to the random samples in input.
        
        Inputs:
            - :``smp``: :
                (``records.Samples`` instance) random samples
        
        .. Note::
            This method assign attribute ``smp`` to the ``mf3.MT`` section
        """
        from sandy.records import Samples
        nsmp = smp.shape[1]
        self.smp = self.xs.init_smp(nsmp, Samples)
        pert = smp.interp(self.xs.xx) + 1
        self.smp *= pert
        self.smp.replace_negative()
    
    def write(self, ismp=None):
        r"""
        Return the data in this ``MT`` section as a list of strings according 
        to the ``ENDF-6`` format.
        
        Inputs:
            - :``ismp``: :
                (integer) sample number
        
        Outputs:
            - :``text``: :
                (list of strings) ``MT`` section data in ``ENDF-6`` format
        """
        from sandy.endf import ENDF, za, awr
        text = []
        text += ENDF.write_cont(za, awr, 0, self.l2, 0, 0)
        if ismp is not None and hasattr(self, 'smp'):
            text += self.smp.write_endf(c1=self.qm, c2=self.qi, l2=self.lr, ismp=ismp)
        else:
            text += self.xs.write_endf(c1=self.qm, c2=self.qi, l2=self.lr)
        text = ENDF.add_control(text, self.mat, self.mf, self.mt)
        return text

    def fig_samples(self, filename=None):
        """
        Plot samples standard deviation against ``MF33`` standard deviation 
        of the first y-axis.
        Plot the cross section on the second y-axis.
        """
        import matplotlib.pyplot as plt
        from os.path import join
        from sandy.plot import save
        from sandy.sandy_input import outdir
        if not hasattr(self, 'smp'):
            return
        fig, ax = plt.subplots()
        if hasattr(self, 'cov'):
            std = self.cov.Std
            std *= 100
            std.plot(ax, step=True, color='r', label='ENDF-6')
        rstd = self.smp.Rstd
        rstd *= 100
        rstd.plot(ax, step=True, color='b', label='samples')
        ax.set_title(self.prefix + "{} samples".format(self.smp.nsmp))
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.set_xlabel('energy (eV)')
        ax.set_ylabel('stdev (%)')
        mask = self.smp.Mean > 0
        xmin, xmax = self.smp.xx[mask][0], self.smp.xx[-1]
        ax.set_xlim([xmin, xmax])
        ymin, ymax = ax.get_ylim()
        ax.set_ylim([0, min(ymax,100)])
        ax.grid()
        ax.legend(loc='best', fancybox=True, shadow=True)
        ax2 = ax.twinx()
        self.xs.plot(ax2, alpha=0.5, color='0.75')
        ax2.set_ylabel(r'cross section (b)', color='0.75')
        ax2.tick_params('y', colors='0.75')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        if filename:
            directory = join(outdir, "plots-" + filename)
            figname = join(directory, 'samples_mf{}_mt{}'.format(self.mf, self.mt))
            save(figname, fig, 'pdf')
        else:
            fig.show()



class MTRedundant(MT):
    """
    ``MT`` section of ``MF3`` for redundant cross sections.
    """
    
    @property
    def daughters(self):
        r"""
        Set of daughters.
        Daughters correspond to the intersection of the redundant cross 
        section's components and the reactions actually present in the 
        ``mf3.MF`` object.
        """
        return self._daughters
    
    @daughters.setter
    def daughters(self, keys):
        """
        If in the set of daughters exist some redundant cross section, the 
        we make sure to remove from the set of daughters also the components 
        of the redundant cross sections that were found.
        
        .. Note:
            Assume that the set of daughters of ``MT1`` is::
                [2, 4, 16, 17, 51, 52, 53, 91, 102]

            ``MT4`` is present in the set of daughters but it is itself a 
            redundant cross section with daughters::
                [51, 52, 53, 91]
            
            The set of daughters for ``MT1`` must be corrected and the 
            daughters of ``MT4`` must be removed, as ::
                [2, 4, 16, 17, 102]
        """
#        import pytest
#        pytest.set_trace()
        daughters = set(keys) & self.__class__.components
        to_remove = set()
        for mt in sorted(daughters & redundant_xs, reverse=True):
            to_remove |= eval("MT{}.components".format(mt))
        self._daughters = daughters - to_remove
    

    def SUM(self, *items):
        r"""
        Reconstruct *redundant* cross section and samples as the sum of its 
        components.
        
        Inputs:
            - items :
                (list) list of ``mf3.MT`` instances
        
        .. Note::
            This method assign attribute ``smp`` to the ``mf3.MT`` section
        """
        from sandy.records import Tab1, Samples
        mts = sorted([item.mt for item in items])
        logging.debug(self.prefix + "Reconstruct reaction as the sum of MT{}".format(mts))
        self.xs = Tab1.SUM([ item.xs for item in items ])
        # Do the same for samples
        if np.any([ hasattr(item, "smp") for item in items ]):
            logging.debug(self.prefix + "Reconstruct samples as the sum of MT{}".format(mts))
            smp = Samples.SUM([ item.smp if hasattr(item, "smp") else item.xs for item in items ])
            self.smp = smp



class MT1(MTRedundant):
    """
    ``MT1`` section of ``MF3`` (total cross section).
    """
    components = set(range(2,850))



class MT3(MTRedundant):
    """
    ``MT3`` section of ``MF3`` (non-elastic cross section).
    """
    components = set(range(4,850))



class MT107(MTRedundant):
    r"""
    ``MT107`` section of ``MF3`` (:math:`(n,\\alpha)` cross section).
    """
    components = set(range(800,850))



class MT106(MTRedundant):
    r"""
    ``MT106`` section of ``MF3`` (:math:`(n,He^{3})` cross section).
    """
    components = set(range(750,800))



class MT105(MTRedundant):
    r"""
    ``MT105`` section of ``MF3`` (:math:`(n,t)` cross section).
    """
    components = set(range(700,750))



class MT104(MTRedundant):
    r"""
    ``MT104`` section of ``MF3`` (:math:`(n,d)` cross section).
    """
    components = set(range(650,700))



class MT103(MTRedundant):
    r"""
    ``MT103`` section of `MF3` (:math:`(n,p)` cross section).
    """
    components = set(range(600,650))



class MT101(MTRedundant):
    """
    ``MT101`` section of `MF3` (neutron disappearance cross section).
    """
    components = set(range(102,118))



class MT18(MTRedundant):
    """
    ``MT18`` section of ``MF3`` (fission cross section).
    """
    components = {19,20,21,38}



class MT27(MTRedundant):
    """
    ``MT27`` section of ``MF3`` (absorption cross section).
    """
    components = {18, 101}



class MT4(MTRedundant):
    """
    ``MT4`` section of ``MF3`` (total inelastic scattering cross section).
    """
    components = set(range(51,92))