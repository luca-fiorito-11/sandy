#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 08:26:25 2017

@author: lfiorito
"""

import sys
import logging



class Perturbation():
    """
    Class dedicated to the perturbation procedure.
    Tabulated cross-section are multiplied by a perturbation coefficient.
    """
    
    def __repr__(self):
        return "<PERTURBATION>"
    
    def __init__(self):
        import sandy.sandy_input
        sandy.sandy_input.options['mf'] = [3]
        self.F = self.get_endf(self.endf)
    
    @property
    def mt(self):
        """
        Input ``mt`` keyword. It must be a list.
        """
        return self.get_keyword('mt', err=True)
    
    @property
    def endf(self):
        """
        Name of the input ``ENDF-6`` file.
        """
        return self.get_keyword('endf', err=True)
    
    @property
    def pert(self):
        """
        List of perturbation coefficients
        """
        return self.get_keyword('pert', err=True)
    
    def get_endf(self, filename):
        from sandy.endf import ENDF
        file = ENDF(filename)
        lrp = file.get_line(1).CONT.l1
        if lrp != 2:
            # Convert to PENDF
            return file.read_pendf(file)
        else:
            file.process()
            return file

    def get_keyword(self, key, err=False):
        from sandy.sandy_input import options
        if key in options:
            return options[key]
        elif err:
            logging.error("ERROR : missing input keyword '{}'".format(key))
            sys.exit()
    
    def run(self):
        r"""
        Run the *perturbation* algorithm.
        For each perturbation coefficient in the input list, perturb the cross 
        sections and write the perturbed output file.
        """
        from sandy.sandy_input import outdir
        from os.path import join
        for pert in self.pert:
            outname = self.F.name + '.perturb_{}'.format(pert)
            for mt in sorted(self.mt):
                if mt not in self.F[3]:
                    logging.warn("PERTURBATION : MT{} is missing in the ENDF-6 files".format(mt))
                    continue
                self.F[3][mt].xs += self.F[3][mt].xs * pert
            self.F[3].sum_rules()
            outfile = join(outdir, outname)
            logging.info("PERTURBATION : Writing file '{}'".format(outfile))
            self.F.write_to_file(outfile)



def preprocessing():
    from sandy.sandy_input import process_input
    process_input()
    logging.info("PERTURBATION : Run perturbation module")
    Perturbation()
    logging.info("CALCULATION COMPLETED")



if __name__ == "__main__":
    preprocessing()
