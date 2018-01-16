# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:16:43 2017

@author: lfiorito
"""
import sys
import logging
from sandy.records import Tab1


class MT(dict):
    
    def __init__(self, DICT):
        if DICT is not None:
            super().__init__(DICT)
            
    @property
    def beg(self):
        if 'beg' in self:
            return float(self['beg'])

    @property
    def end(self):
        if 'end' in self:
            return float(self['end'])


class Replace():
    
    @property
    def prefix(self):
        return "REPLACE : "
        
    @property
    def mt(self):
        DICT = self.get_keyword('replace', err=True)
        if hasattr(DICT, "items"):
            return { k : MT(v) for k,v in DICT.items() }
        else:
            return { DICT : MT({}) }
    
    def get_endf(self, key):
        from sandy.endf import ENDF
        filename = self.get_keyword(key, err=True)
        endf = ENDF(filename)
        endf.process()
        lrp = endf.get_line(1).CONT.l1
        if lrp != 2:
            # Convert to PENDF
            pendf = ENDF.read_pendf(endf)
        else:
            pendf = endf
        return endf, pendf

    def get_keyword(self, key, err=False):
        from sandy.sandy_input import options
        if key in options:
            return options[key]
        elif err:
            logging.error("ERROR : missing input keyword '{}'".format(key))
            sys.exit()
    
    def run(self):
        r"""
        """
        from sandy.sandy_input import outdir, plot
        from os.path import join
        logging.info("\n"+80*"=")
        logging.info("{:^80}".format("REPLACE"))
        logging.info(80*"=")
        From, FROM = self.get_endf('from')
        To, TO = self.get_endf('to')
        outname = TO.name + '.replace'
        outfile = join(outdir, outname)
        logging.debug(80*"+")
        logging.debug("{:^80}".format(" Replace cross sections "))
        logging.debug(80*"+")
        for mt,MT in sorted(self.mt.items()):
            if mt not in FROM[3] or mt not in TO[3]:
                logging.warn(self.prefix + "MT{} is missing in one or both ENDF-6 files".format(mt))
                continue
            xs_from = FROM[3][mt].xs
            xs_to = TO[3][mt].xs.copy()
            xs = xs_from.reinit(xs_to.xx)
            mask = xs.get_section(MT.beg, MT.end)
            if MT.beg is not None:
                beg = MT.beg
            else:
                beg = xs.xx[mask][0]
                logging.info(TO[3][mt].prefix + "'beg' not given. Set to {:.3} eV".format(beg))
            if MT.end is not None:
                end = MT.end
            else:
                end = xs.xx[mask][-1]
                logging.info(TO[3][mt].prefix + "'end' not given. Set to {:.3} eV".format(end))
            msg = r"Replace data for {:.3} <= E <= {:.3} eV"
            logging.info(TO[3][mt].prefix + msg.format(beg, end))
            TO[3][mt].xs[mask] = xs[mask]
            if plot:
                # Plotting
                labels = ['FROM', 'TO', 'REPLACE']
                title = r"MT{} : replace data for ${:.3} \leq E \leq {:.3} eV$".format(mt, beg, end)
                name = outname + '_mt{}'.format(mt)
                Tab1.figs(xs_from, xs_to, TO[3][mt].xs, 
                          xlabel='energy (eV)', ylabel='cross section (b)', 
                          labels=labels, title=title, name=name)
        TO[3].sum_rules()
        if 152 in TO[2]:
            logging.debug(self.prefix + "Copy section MF2/MT152 from PENDF")
            To[2].update({ 152 : TO[2][152]})
        logging.debug(self.prefix + "Copy section MF3 from PENDF")
        To.update({ 3 : TO[3]})        
        logging.info("REPLACE : Writing file '{}'".format(outfile))
        To.write_to_file(outfile)



def preprocessing():
    from sandy.sandy_input import process_input
    options = process_input()
    logging.info("REPLACE : Run replace module")
    Replace(**options)
    logging.info("CALCULATION COMPLETED")



if __name__ == "__main__":
    preprocessing()
