import logging
import sys

from sandy.formats import *
from sandy.settings import *
from sandy.decay import *
from sandy.njoy import *
from .core import *
import sandy.tools

class ShutdownHandler(logging.Handler):
    """Trigger exit on errors.
    """
    def emit(self, record):
        logging.shutdown()
        sys.exit(1)

class DuplicateFilter(object):
    """Define a filter which keeps track of what was logged, and attach it to 
    your logger for the duration of a loop.
    """
    def __init__(self):
        self.msgs = set()

    def filter(self, record):
        rv = record.msg not in self.msgs
        self.msgs.add(record.msg)
        return rv

class Error(Exception):
    pass
    
FORMAT = '%(levelname)s:  %(message)s'
logging.basicConfig(format=FORMAT)
logging.getLogger().setLevel(logging.INFO)
logging.getLogger().addHandler(ShutdownHandler(level=40))
#logging.getLogger().addFilter(DuplicateFilter())

__version__ = '0.9.0'
