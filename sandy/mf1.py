# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 13:35:32 2017

@author: lfiorito
"""

import logging
import sys
from sandy.records import List
import numpy as np
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard


class PolyNubar(List):
    """
    Nubar (:math:`\\overline{\\nu}`) given in polynomial form.
    Values of :math:`\\overline{\\nu}` are tabulated in the form of a list 
    of coefficients provided for the following polynomial expansion

    .. math::
        
        \\overline{\\nu}(E) = \\sum_{n=1}^{NC} C_n E^{n-1}
    
    where:
        - :math:`\\overline{\\nu}(E)` :
            is the average number of neutrons per fission produced by neutrons 
            of incident energy E (eV)
        - :math:`C_n` :
            are the *n-th* coefficient
        - :math:`NC` :
            is the length of the list and corresponds to the number of terms 
            in the polynomial
    """
    
    def convert_to_tab(self, emin=1e-5, emax=2e7, npoints=50):
        """
        Convert the nubar function from polynomial expansion to tabulated.
        
        Outputs:
            - nubar :
                (`Tab1` instance) tabulated nubar
        """
        from sandy.records import Tab1
        MIN = np.log10(emin)
        MAX = np.log10(emax)
        xx = np.logspace(MIN, MAX, npoints)
        yy = np.sum([ c*xx**i for i,c in enumerate(self) ], axis=0)
        nubar = Tab1(xx, yy, intscheme=[2])
        return nubar



class MF(MFstandard):
    
    def __init__(self, endf):
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Processing ENDF-6 section MF1 : GENERAL INFORMATION "))
        logging.debug(80*"+")
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            if mt == 452:
                self[mt] = MT452(endf)
            elif mt == 455 or mt == 456:
                self[mt] = MT(endf)
            elif mt != 451:
                self[mt] = MTstandard(endf)
            endf.next() # send
        endf.next() # fend
        self._set_daughters()
    
    @property
    def union_grid(self):
        r"""
        Union of the energy grids for all the neutron multiplicity sections.
        """
        from sandy.functions import union_grid
        grids = [mtsec.nubar.xx for mt,mtsec in self.items() if mt in [452, 455, 456]]
        ug = union_grid(*grids)
        return ug

    def _set_daughters(self):
        """
        Assign *daughters* to the total nubar.
        For further info, see docstring of attribute 
        ``mf3.MTRedundant.daughters``.
        
        .. Important:
            This method assigns attribute ``daughters`` to the ``mf1.MT`` 
            objects
        """
        for mt in self:
            if hasattr(self[mt], "components"):
                self[mt].daughters = self.keys()

    def MT_NUBAR(self):
        r"""
        Generator that yields nubar sections (``MT452``, ``MT455``, ``MT456``).
        """
        mts = {452, 455, 456} & self.keys()
        DICT = {mt : self[mt] for mt in mts}
        return DICT.items()

    def is_unionized(self, err=True):
        """
        Check if all the nubar are *unionized* (defined on a common union grid).
        
        Inputs:
            - :``err``: :
                (boolean) flag that raises error or not
        
        Outputs:
            - :``is_unionized``: : 
                (boolean) ``True`` if the cross sections are *unionized*, 
                otherwise ``False``
        """
        is_unionized = []
        for mt,mtsec in sorted(self.MT_NUBAR()):
            check = mtsec.nubar.xx.tolist() == self.union_grid.tolist()
            if not check:
                raise NotImplementedError("nubar values cannot be summed because 'MT{}' is not unionized".format(mt))
            is_unionized.append(check)
        is_unionized = np.all(is_unionized)
        return is_unionized

    def sum_rules(self):
        r"""
        Reconstruct total nubar ``MT452`` as the sum of ``MT455`` and 
        ``MT456``.
        This method works only if all the tabulated functions are unionized.
        """
        self.is_unionized()
        logging.info(self.prefix + "Reconstruct redundant neutron multiplicities")
        self[452].SUM(self[455], self[456]) 

    def unionize(self, xx=[]):
        r"""
        Interpolate all the neutron multiplicity tabulated functions over 
        a common grid.
        
        Inputs:
            - :``xx``: :
                (iterable) energy grid, by default take the union grid of the 
                existing nubar sections
        """
        from sandy.records import Grid
        logging.debug(self.prefix + "Reconstruct nubar tabulated functions over a common energy grid")
        for mt,mtsec in self.MT_NUBAR():
            grid = self.union_grid if len(xx) < 2 else Grid(xx)
            mtsec.nubar = mtsec.nubar.reinit(grid)

    def add_samples(self, mf31):
        """
        Add samples from ``MF31`` to ``MF1``.
        
        Inputs:
            - :``mf31``: :
                (``mf31.MF`` instance) dictionary containing ``MF31``
        
        .. Note::
            This method is analogous to ``mf3.MF.add_samples``, except that 
            the nubar functions are interpolated over a union grid.
        """
        from sandy.functions import union_grid
        extra_points_mf31 = mf31.union_grid.add_extra_points()
        ug = union_grid(self.union_grid, extra_points_mf31)
        self.unionize(ug)
        for mt, mtsec in sorted(mf31.items(), reverse=True):
            if mt not in self:
                logging.warn(mtsec.prefix + "Section MF{} MT{} was not found".format(self.mf, mt))
                continue
            smp = mtsec.smp
            cov = mtsec[mt]
            logging.debug(mtsec.prefix + "Assign samples to MF{} MT{}".format(self.mf, mt))
            self[mt].nubar = self[mt].nubar.add_points(ug, prefix="MF{} MT{} : ".format(self.mf, mt))
            self[mt]._add_samples(smp)
            self[mt].cov = cov
            if hasattr(self[mt], "daughters"):
                # Check that no daughter has samples
                if not np.any([ hasattr(self[x], "smp") for x in self[mt].daughters ]):
                    for dau in self[mt].daughters:
                        logging.debug(mtsec.prefix + "Assign samples to MF{} MT{}".format(self.mf, dau))
                        self[dau].nubar = self[dau].nubar.add_points(ug, prefix="MF{} MT{} : ".format(self.mf, dau))
                        self[dau]._add_samples(smp)
                        self[dau].cov = cov
        self.sum_rules()

    def plot_samples(self, filename):
        """
        Call the samples-plotting method of each ``MT`` section containing 
        nubar.
        
        Inputs:
            - :``filename``: :
                (string) name of the ``ENDF-6`` file
        """
        for mt,mtsec in sorted(self.MT_NUBAR()):
            mtsec.fig_samples(filename)



class MT(MTstandard):
    """
    ``MT452``, ``MT455`` or ``MT456`` section of ``MF1``.
    These sections are used to store tabulated or polynimal representations 
    of total, delayed and prompt fission neutron multiplicities, respectively.
    """
    
    def __init__(self, endf):
        from sandy.records import Tab1
        self.mat, self.mf, self.mt, ns = endf.read_control()
        za, awr, self.ldg, self.lnu = endf.read_cont([4,5])
        if self.mt == 455:
            self.lambdas = List.read_endf(endf)
        if self.lnu == 1: # (endf/b-vii.1 Ra223)
            nubar = PolyNubar.read_endf(endf)
            self.nubar = nubar.convert_to_tab(emax=endf.INFO.emax)
            self.lnu = 2
        elif self.lnu == 2:  # (jendl-4.0 U235)
            self.nubar = Tab1.read_endf(endf)
    
    @property
    def ldg(self):
        r"""
        Flag indicating energy dependence of delayed-group constants.
        - :``ldg=0``: :
            decay constants are energy-independent
        - :``ldg=1``: :
            decay constants are energy-dependent
        
        ..Important::
            In ``ENDF/B-VII.1``, ``JEFF-3.2``, ``JENDL-4.0`` and 
            ``TENDL-2015`` value ``ldg=1`` was never found.
            Hence, this capability was not included in ``SANDY``.
        """
        return self._ldg

    @ldg.setter
    def ldg(self, ldg):
        if ldg != 0 and self.mt != 455:
            logging.error(self.prefix + "Attribute 'ldg!=0' is acceptable only for MT455, not for MT{}}".format(self.mt))
            sys.exit()
        elif ldg == 1 and self.mt == 455:
            logging.error(self.prefix + "SANDY cannot yet process energy-dependent delayed-group constants")
            sys.exit()
        self._ldg = ldg
    
    @property
    def lnu(self):
        r"""
        Indicates what representation of nubar has been used.
        - :``lnu=1``: :
            polynomial representation
        - :``lnu=2``: :
            tabulated representation

        ..Important::
            In ``ENDF/B-VII.1``, ``JEFF-3.2``, ``JENDL-4.0`` and 
            ``TENDL-2015`` value ``lnu=2`` was never found for prompt or 
            delayed nubar.
            Hence, this capability was not included in ``SANDY``.
        """
        return self._lnu

    @lnu.setter
    def lnu(self, lnu):
        if lnu == 1 and self.mt != 452:
            logging.error(self.prefix + "SANDY cannot yet process polynomial representation")
            sys.exit()
        self._lnu = lnu

    def _add_samples(self, smp):
        r"""
        Perturb nubar according to the random samples in input.
        
        Inputs:
            - :``smp``: :
                (``records.Samples`` instance) random samples
        
        .. Note::
            This method assign attribute ``smp`` to the ``mf1.MT`` section
        
        .. Note::
            This method is analogous to ``mf3.MT._add_samples``
        """
        from sandy.records import Samples
        nsmp = smp.shape[1]
        self.smp = self.nubar.init_smp(nsmp, Samples)
        pert = smp.interp(self.nubar.xx) + 1
        self.smp *= pert
        self.smp.replace_negative()

    def write(self, ismp=None):
        r"""
        Return the data in this ``MT`` section as a list of strings according 
        to the ``ENDF-6`` format.
        
        Inputs:
            - :``ismp``: :
                (integer) sample number, default is ``None``
        
        Outputs:
            - :``text``: :
                (list of strings) ``MT`` section data in ``ENDF-6`` format
        """
        from sandy.endf import ENDF, za, awr
        text = []
        text += ENDF.write_cont(za, awr, self.ldg, self.lnu, 0, 0)
        if self.mt == 455:
            text += self.lambdas.write_endf()
        if self.lnu == 1 and self.mt == 452:
            text += self['c'].write_endf()
        elif self.lnu == 2:  # (jendl-4.0 U235)
            if ismp is not None and hasattr(self, 'smp'):
                text += self.smp.write_endf(ismp=ismp)
            else:
                text += self.nubar.write_endf()
        text = ENDF.add_control(text, self.mat, self.mf, self.mt)
        return text
    
    def fig_samples(self, filename=None):
        """
        Plot samples standard deviation against ``MF31`` standard deviation 
        of the first y-axis.
        Plot the neutron multiplicity on the second y-axis.
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
        self.nubar.plot(ax2, alpha=0.5, color='0.75')
        ax2.set_ylabel(r'nubar (-)', color='0.75')
        ax2.tick_params('y', colors='0.75')
        ax2.set_xscale('log')
        ax2.set_yscale('linear')
        if filename:
            directory = join(outdir, "plots-" + filename)
            figname = join(directory, 'samples_mf{}_mt{}'.format(self.mf, self.mt))
            save(figname, fig, 'pdf')
        else:
            fig.show()



class MT452(MT):
    """
    ``MT452`` section of ``MF1`` (total nubar).
    """
    components = {455, 456}
    
    @property
    def daughters(self):
        r"""
        See docstring of ``mf3.MTRedundant.daughters``.
        """
        return self._daughters
    
    @daughters.setter
    def daughters(self, keys):
        self._daughters = set(keys) & self.__class__.components
    
    def SUM(self, *items):
        r"""
        Reconstruct *redundant* cross section and samples as the sum of its 
        components.
        
        Inputs:
            - items :
                (list) list of ``mf3.MT`` instances
        
        .. Note::
            This method assign attribute ``smp`` to the ``mf3.MT`` section

        .. Note::
            This method is analogous to ``mf3.MTRedundant.SUM``
        """
        from sandy.records import Tab1, Samples
        mts = sorted([item.mt for item in items])
        logging.debug(self.prefix + "Reconstruct reaction as the sum of MT{}".format(mts))
        self.nubar = Tab1.SUM([ item.nubar for item in items ])
        # Do the same for samples
        if np.any([ hasattr(item, "smp") for item in items ]):
            logging.debug(self.prefix + "Reconstruct samples as the sum of MT{}".format(mts))
            smp = Samples.SUM([ item.smp if hasattr(item, "smp") else item.xs for item in items ])
            self.smp = smp



class Info:
    """
    ``MF1`` ``MT451`` info section
    """
    
    def __init__(self, endf):
        import sandy.endf
        self.mat, self.mf, self.mt, ns = endf.read_control()
        sandy.endf.za, sandy.endf.awr, self.lrp, self.lfi, self.nlib, self.nmod = endf.read_cont()
        self.elis, self.sta, self.lis, self.liso, self.nfor = endf.read_cont([4])
        self.awi, self.emax, self.lrel, self.nsub, self.nver = endf.read_cont([3])
        self.temp, self.ldrv, nwd, nxc = endf.read_cont([1,3])
        self.info = [ endf.read_text() for i in range(nwd) ]
        self.records = Records.read_records(nxc, endf)

    @property
    def izam(self):
        """
        .. math:
            ZZZ + AAA + I
        
        where:
            - :math:`ZZZ` is the atomic number in 3 digits
            - :math:`AAA` is the mass number in 3 digits
            - :math:`I` is the isomeric state in 1 digit
        """
        from sandy.endf import za
        return int(za*10) + self.liso
    
    @property
    def lrp(self):
        r"""
        Flag indicating whether resolved and/or unresolved resonance 
        parameters are given in `MF2`.
        
        - lrp=-1 :
            no File 2 is given (not allowed for incident neutrons)
        - lrp=0 :
            no resonance parameter data are given, but a File 2 is present 
            containing the effective scattering radius
        - lrp=1 :
            resolved and/or unresolved parameter data are given in File 2 and
            cross sections computed from them must be added to background 
            cross sections given in File 3
        - lrp=2 :
            parameters are given in File 2, but cross sections derived from 
            them are not to be added to the cross sections in File 3.
            This option is to be used for derived files only and is typical 
            in the so-called `PENDF` files, in which the cross sections are 
            already reconstructed from the resonances parameters and written 
            in File 3
        """
        return self._lrp
    
    @lrp.setter
    def lrp(self, lrp):
        self._lrp = lrp

    @property
    def lfi(self):
        """
        Flag indicating whether this material fissions.
        
        - lfi=0 :
            this material does not fission
        - lfi=1 :
            this material fissions
        """
        return self._lfi
    
    @lfi.setter
    def lfi(self, lfi):
        self._lfi = lfi

    @property
    def nlib(self):
        r"""
        Library identifier (e.g. `nlib=0` for `ENDF/B`). 
        Additional values have been assigned to identify other libraries using 
        `ENDF-6` format.
        """
        return self._nlib
    
    @nlib.setter
    def nlib(self, nlib):
        self._nlib = nlib

    @property
    def nmod(self):
        r"""
        Modification number for this material.
        - nmod=0 :
            evaluation converted from a previous version
        - nmod=1 :
            new or revised evaluation for the current library version
        - nmod>1 :
            for successive modifications
        """
        return self._nmod
    
    @nmod.setter
    def nmod(self, nmod):
        self._nmod = nmod

    @property
    def elis(self):
        r"""
        Excitation energy of the target nucleus relative to 0.0 for the ground 
        state.
        """
        return self._elis
    
    @elis.setter
    def elis(self, elis):
        self._elis = elis

    @property
    def sta(self):
        r"""
        Target stability flag.
        - sta=0 :
            stable nucleus
        - sta=1 :
            unstable nucleus.
            If the target is unstable, radioactive decay data should be given 
            in the decay data sub-library (nsub=4).
        """
        return self._sta
    
    @sta.setter
    def sta(self, sta):
        self._sta = sta

    @property
    def lis(self):
        r"""
        State number of the target nucleus.
        The ground state is indicated by `lis=0`.
        """
        return self._lis
    
    @lis.setter
    def lis(self, lis):
        self._lis = lis

    @property
    def liso(self):
        r"""
        Isomeric state number.
        The ground state is indicated by `liso=0`.
        `lis` is greater than or equal to `liso`.
        """
        return self._liso
    
    @liso.setter
    def liso(self, liso):
        self._liso = liso

    @property
    def nfor(self):
        r"""
        Library format.
        `nfor=6` for all libraries prepared according to the specifications 
        given in the `ENDF-6` manual.
        """
        return self._nfor
    
    @nfor.setter
    def nfor(self, nfor):
        if nfor != 6:
            logging.error("MF1 MT451 : File is not written in `ENDF-6` format")
            sys.exit()
        self._nfor = nfor

    @property
    def awi(self):
        r"""
        Mass of the projectile in neutron mass units
        For incident photons or decay data sub-libraries, use `awi=0`.
        """
        return self._awi
    
    @awi.setter
    def awi(self, awi):
        self._awi = awi

    @property
    def emax(self):
        r"""
        Upper limit of the energy range for evaluation.
        """
        return self._emax
    
    @emax.setter
    def emax(self, emax):
        self._emax = emax

    @property
    def lrel(self):
        r"""
        Library release number; for example, `lrel=2` for the `ENDF/B-VI.2` 
        library.
        """
        return self._lrel
    
    @lrel.setter
    def lrel(self, lrel):
        self._lrel = lrel

    @property
    def nsub(self):
        r"""
        Sub-library number.
        """
        return self._nsub
    
    @nsub.setter
    def nsub(self, nsub):
        self._nsub = nsub

    @property
    def nver(self):
        r"""
        Library version number; for example, `nver=7` for version `ENDF/B-VII`.
        """
        return self._nver
    
    @nver.setter
    def nver(self, nver):
        self._nver = nver

    @property
    def temp(self):
        r"""
        Target temperature (Kelvin) for data that have been generated by 
        Doppler broadening.
        For derived data only; use `temp=0.0` for all primary evaluations.
        """
        return self._temp
    
    @temp.setter
    def temp(self, temp):
        self._temp = temp

    @property
    def ldrv(self):
        r"""
        Special derived material flag that distinguishes between different 
        evaluations with the same material keys (i.e., `mat`, `nmod`, `nsub`).
        - ldrv=0 :
            primary evaluation
        - ldrv>0 :
            special derived evaluation (for example, a dosimetry evaluation
            using sections (`MT`) extracted from the primary evaluation).
        """
        return self._ldrv
    
    @ldrv.setter
    def ldrv(self, ldrv):
        self._ldrv = ldrv

    def write(self, fendf, ismp=None):
        """
        Return MT section's text as a list of strings.
        
        Outputs:
            - text : 
                list of strings containing the MT section
        """
        from sandy.endf import ENDF, za, awr
        text = []
        text += ENDF.write_cont(za, awr, self.lrp, self.lfi, self.nlib, self.nmod)
        text += ENDF.write_cont(self.elis, self.sta, self.lis, self.liso, 0, self.nfor)
        text += ENDF.write_cont(self.awi, self.emax, self.lrel, 0, self.nsub, self.nver)
        infolines = self.write_info()
        self.update_records(fendf)
        recordslines = self.records.write()
        nwd = len(infolines)
        nxc = len(recordslines)
        records = (self.temp, 0, self.ldrv, 0, nwd, nxc)
        text.extend(ENDF.write_cont(*records))
        text.extend(infolines)
        text.extend(recordslines)
        text = ENDF.add_control(text, self.mat, self.mf, self.mt)
        return text
    
    def update_records(self, endf):
        records = []
        for mf,mfsec in sorted(endf.items()):
            for mt,mtsec in sorted(mfsec.items()):
                if (mf,mt) != (1,451):
                    records.append((mf, mt, len(mtsec.text), 0))
        len1451 = 4 + len(self.info) + len(records) + 1
        records = [(1, 451, len1451, 0)] + records
        self.records = Records(*records)
    
    def write_info(self):
        text = [ info for info in self.info ]
        return text



class Records(dict):

    @classmethod
    def read_records(cls, nxc, endf):
        """
        Read number of records in the directory for this material.
        
        Each section (`MT`) in the material has a corresponding line in the 
        directory that contains `MF`, `MT`, `NC`, and `MOD`.
        
        `NC` is a count of the number of records in the section (not 
        including `SEND`), and `MOD` is the modification flag.
        
        Below is an examples of the records for the U238 file from JENDL-3.3::
            
            1        451        554          49237 1451  410
            1        452          6          19237 1451  411
            1        455          7          49237 1451  412
            1        456          5          19237 1451  413
            2        151       2538          39237 1451  414
            3          1        119          19237 1451  415
            3          2        100          49237 1451  416
            3          4         46          49237 1451  417
            3         16         14          49237 1451  418
            3         17         10          49237 1451  419
            3         18         60          49237 1451  420
            3         37          5          19237 1451  421
            3         51         37          19237 1451  422
            3         52         33          19237 1451  423
            3         53         31          19237 1451  424
            3         54         26          49237 1451  425
            3         55         25          49237 1451  426
            3         56         24          49237 1451  427
            3         57         24          19237 1451  428
            3         58         24          19237 1451  429
        """
        records = [ endf.read_cont([0,1]) for i in range(nxc) ]
        inst = cls(*records)
        return inst
        
    def __init__(self, *records):
        """
        Reset records dictionary in `MF1` `MT451` based on input parameters.
        
        Inputs:
            - records :
                series of records as positional arguments.
                Each record is a sequence of four parameters: mf, mt, nc, mod.
                If mod is None, set it to zero.
        """
        for mf, mt, nc, mod in records:
            if mf not in self:
                self[mf] = {}
            self[mf][mt] = nc#{'nc' : nc, 'mod' : mod}

    def write(self):
        record = "{:22}{:>11}{:>11}{:>11}{:>11}"
        text = []
        for mf in sorted(self):
            for mt in sorted(self[mf]):
                nc = self[mf][mt]
                mod = 0
                string = record.format('', mf, mt, nc, mod)
                text.append(string)
        return text