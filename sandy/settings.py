# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:21:37 2018

@author: fiorito_l
"""
def init(ARGS=None):
    import os
    import argparse
    global args
    parser = argparse.ArgumentParser(description='Run SANDY')
    parser.add_argument('endf6',
                        help="ENDF-6 format file")
    parser.add_argument('--samples',
                        type=int,
                        default=100,
                        help="Number of samples")
    parser.add_argument('--outdir',
                        default=os.getcwd(),
                        help="Target directory were outputs are stored. If it does not exist it will be created")
    parser.add_argument('--processes',
                        type=int,
                        default=None,
                        help="Number of worker processes. By default, the number returned by os.cpu_count() is used")
    parser.add_argument('--pendf',
                        help="PENDF format file")
    parser.add_argument('-nw',
                        '--no-write',
                        help="Disable writing")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="SANDY's version")
    args = parser.parse_args(args=ARGS)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)
