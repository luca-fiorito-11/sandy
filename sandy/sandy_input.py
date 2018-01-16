# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:48:11 2017

@author: lfiorito
"""
import sys
import logging
from argparse import ArgumentParser


global plot, write, outdir, options
plot = False
write = True
outdir = None
options = {}
title =  "\n".join((" _______  _______  __    _  ______   __   __ ",
                    "|       ||   _   ||  |  | ||      | |  | |  |",
                    "|  _____||  |_|  ||   |_| ||  _    ||  |_|  |",
                    "| |_____ |       ||       || | |   ||       |",
                    "|_____  ||       ||  _    || |_|   ||_     _|",
                    " _____| ||   _   || | |   ||       |  |   |  ",
                    "|_______||__| |__||_|  |__||______|   |___|  ",
                    "                                             \n"))


class List(list):
    """
    ``list`` object.
    """

    def __init__(self, item, size=None, TYPE=int):
        LIST = item if hasattr(item, "__len__") else [item]
        super().__init__(LIST)
    
    def _check_size(self, size):
        pass

    def _check_type(self, TYPE):
        try:
            LIST = list(map(TYPE, self))
        except:
            raise NotImplementedError("Cannot convert list into type {}".format(TYPE))
        super().__init__(LIST)


class MF(List):
    r"""
    Class containing the arguments of input key *mf*.
    """

    def __init__(self, item):
        super().__init__(item)
        try:
            self._check_type(int)
        except NotImplementedError:
            logging.error("INPUT : 'mf' argument must be a list of integers")
        if 31 in self:
            self.append(1)
        if 32 in self:
            self.append(2)
        if 33 in self:
            self.append(3)
        if 34 in self:
            self.append(4)
        if 35 in self:
            self.append(5)
        super().__init__(list(set(self)))



class MT(List):
    r"""
    Class containing the arguments of input key *mt*.
    """

    def __init__(self, item):
        super().__init__(item)
        try:
            self._check_type(int)
        except NotImplementedError:
            logging.error("INPUT : 'mt' argument must be a list of integers")



class MODE(str):
    r"""
    Class containing the argument of input key *mode*.
    """
    
    modes = ["sampling", "replace", "perturbation"]
    
    def __init__(self, item):
        try:
            mode = item.lower()
        except:
            logging.error("INPUT : 'mode' argument '{}' is not a string".format(item))
            sys.exit()
        if mode not in self.modes:
            logging.error("INPUT : 'mode' argument '{}' is not a valid mode".format(item))
            sys.exit()



class Options(dict):
    """
    Dictionary of ``SANDY`` input options.
    """
    
    def __init__(self, **kwargs):
        super().__init__()
        for k,v in kwargs.items():
            self[k] = v
        if "mode" not in self:
            logging.error("INPUT : Missing keyword 'mode' in input file")
            sys.exit()
        if self["mode"] == "sampling":
            self._set_mf_sampling()
        elif self["mode"] == "replace":
            self['mf'] = 3
#        self.check_list('mt')
    
    def __setitem__(self, key, item):
        r"""
        Before setting options, convert keys into lowercase.
        Then, create proper instance for each option argument.
        """
        lkey = key.lower()
        if lkey == "mode":
            super().__setitem__(lkey, MODE(item))
        elif lkey == 'mf':
            super().__setitem__(lkey, MF(item))
        elif lkey == 'mt':
            super().__setitem__(lkey, MT(item))
        else:
            super().__setitem__(lkey, item)
    
    def check_list(self, key, TYPE=int):
        r"""
        If option is given, check if its argument is a list.
        If the argument is a scalar, make it a list.
        Optionally, check the type of the list elements.
        """
        if key in self:
            if not hasattr(self[key], "__len__"):
                self[key] = [self[key]]
            if not isinstance(self[key], list):
                logging.error("INPUT : Option '{}' must be a list".format(key))
                sys.exit()
            if TYPE is int:
                for i,item in enumerate(self[key]):
                    try:
                        self[key][i] = int(item)
                    except ValueError:
                        logging.error("INPUT : '{}' component '{}' must be an integer".format(key, item))
                        sys.exit()

    def _set_mf_sampling(self):
        r"""
        Predefine the non-sorted list of ``MF`` values in *mf* option.
        If a covariance ``MF`` is given, add also its corresponding data, 
        e.g if ``MF33`` is present, then add ``MF3``.
        
        The standard *mf* list (provided when *mf* is not set) is:
            * MF1 :
                fission neutron multiplicities
            * MF31 :
                fission neutron multiplicities covariances
            * MF2 :
                resonance parameters
            * MF32 :
                resonance parameters covariances
            * MF3 :
                cross sections
            * MF33 :
                cross sections covariances
            * MF4 :
                angular distributions
            * MF34 :
                angular distributions covariances
            * MF5 :
                energy distributions
            * MF35 :
                energy distributions covariances
        """
        if 'mf' not in self:
            self['mf'] = [1,31,2,32,3,33,4,34,5,35]
        self.check_list('mf')
        for i in (1,2,3,4,5):
            if 30+i in self['mf'] and i not in self['mf']:
                self['mf'].append(i)



class CmdLine(ArgumentParser):
    
    def __init__(self, args=None):
        """
        Parse the command line when running SANDY.
        
        ``usage: sandy_input.py [-h] [-o OUTPUT] [--version] input``
        """
        from os import getcwd
        super().__init__(description='Run SANDY')
        self.add_argument('input',
                          help="SANDY's inputfile")
        self.add_argument('-o', 
                          '--output',
                          default=getcwd(), 
                          help="SANDY's output folder")
        self.add_argument('-l', 
                          '--logging',
                          default="debug",
                          help="set logging level (debug, info, warn)")
        self.add_argument("-p",
                          "--plot",
                          action="store_true",
                          help="enable plotting")
        self.add_argument("-q",
                          "--quiet",
                          action="store_true",
                          help="run in quiet mode")
        self.add_argument("--no-write",
                          action="store_false",
                          help="disable writing")
        self.add_argument('--version',
                          action='version',
                          version='%(prog)s 1.0',
                          help="SANDY's version")
        if args is not None:
            args = args.split()
        self.parse_args(args=args, namespace=self)



def lowercase(din):
    dout = {}
    for k,v in din.items():
        try:
            key = k.lower()
        except:
            logging.error("INPUT : non-recognized keyword '{}'".format(k))
            sys.exit()
        dout[key] = v
    return dout



def get_endf_files():
    from os.path import isdir, isfile, join, abspath, expanduser
    from os import listdir
    if 'endf' not in options:
        logging.error("INPUT : Missing keyword 'endf'")
        sys.exit()
    path = abspath(expanduser(options.pop('endf')))
    if isdir(path):
        files = [ join(path,name.lstrip()) for name in listdir(path) ]
    elif isfile(path):
        files = [path]
    else:
        logging.error(" INPUT : Cannot recognize 'endf' input '{}'. It must be a file or a directory".format(path))
        sys.exit()
    return files



def setup_logging(log="info"):
    r"""
    Setup the logging level, default is `INFO`.
    If requested, set quiet mode by using level 100.
    """
    logs = {"warn"  : logging.WARN,
            "debug" : logging.DEBUG,
            "info"  : logging.INFO,
            "quiet" : 100}
    level = logs[log] if log in logs else logging.INFO
    FORMAT = '%(message)s'
    logging.basicConfig(format=FORMAT, stream=sys.stderr, level=level)



def process_input():
    r"""
    Process ``SANDY``'s input.
    
    Outputs:
        - :``inp``: :
            (dictionary) ``SANDY``'s input options
    """
    from sandy.files import File
    from os.path import exists
    from os import makedirs
    global plot, write, outdir, options
    sys.stdout.write(title)
    if CmdLine().quiet:
        log = "quiet"
    else:
        log = CmdLine().logging.lower()
    setup_logging(log)
    filename = CmdLine().input
    options = Options(**File(filename).load_yaml())
    # Create output directory
    outdir = CmdLine().output
    if not exists(outdir):
        makedirs(outdir)
    plot = CmdLine().plot
    write = CmdLine().no_write