# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 11:10:07 2018

@author: Fiorito_L
"""

import argparse
import os
import json
import logging
import pdb

from ..settings import SandyError
from ..formats import read_formatted_file
from ..utils import is_valid_file

__author__ = "Luca Fiorito"
__all__ = ["replace"]

def _parse(iargs=None):
    """Parse command line arguments for replace option.
    """
    parser = argparse.ArgumentParser(prog="python -m sandy.replace",
#                                     description='Run sampling',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('target',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 or PENDF format file.")
    parser.add_argument('source',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 or PENDF format file.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--sections','-S',
                       type=str,
                       help="sections to copy from <source> to <target>.")
    group.add_argument('--ptable',
                       default=True,
                       action="store_true",
                       help="Add probabiliy tables from MF2/MT153.")
    parser.add_argument('--output', '-O',
                        metavar="OUT",
                        type=str,
                        help="output file (default is '<source>.replace').")
#    parser.add_argument('--reconstruct',
#                        metavar="OUT",
#                        default=True,
#                        action="store_false",
#                        help="reconstruct redundant quantities."),
    return parser.parse_known_args(args=iargs)[0]

def replace(iargs=None):
    """Replace/add MF/MT sections from one ENDF-6 file to another.
    """
    init = _parse(iargs)
    if init.ptable:
        sections = {2 : [152]}
    else:
        sections = json.loads(init.sections)
    target = read_formatted_file(init.target)
    output = target.add_sections(init.source, sections)
    string = output.update_info().write_string()
    outname = os.path.split(init.target)[1] + ".replace" if init.output is None else init.output
    with open(outname, 'w') as f: f.write(string)
    print("output written in '{}'".format(outname))
    
def run():
    try:
        replace()
    except SandyError as exc:
        logging.error(exc.args[0])