#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:03:32 2017

@author: lfiorito
"""

import logging

def preprocessing():
    r"""
    Preprocess ``SANDY``'s inputs and select running mode.
    """
    from sandy.sandy_input import process_input, get_endf_files
    from sandy.replace import Replace
    from sandy.sampling import Sampling
    from sandy.perturbation import Perturbation
    process_input()
    from sandy.sandy_input import options
    if options['mode'] == 'replace':
        Replace().run()
    elif options['mode'] == 'perturbation':
        Perturbation().run()
    elif options['mode'] == 'sampling':
        files = get_endf_files()
        for file in files:
            Sampling(file)
    logging.info("CALCULATION COMPLETED")

if __name__ == "__main__":
    preprocessing()
