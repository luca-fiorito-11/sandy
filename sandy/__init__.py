import logging
import sys

from sandy.formats import read_formatted_file
from sandy.settings import *
from sandy.constants import *
from sandy.decay import *
from sandy.energy_grids import *
from sandy.errorr import *
from sandy.groupr import *
from sandy.fy import *
from sandy.h5file import *
from sandy.pert import *
from sandy.pfns import *
from sandy.processing import *
from sandy.zam import *
from sandy.njoy import *
from sandy.sections import *
from sandy.shared import *
from sandy.utils import *
from sandy.core import *
from sandy.sampling import *
import sandy.mcnp
import sandy.tools
import sandy.shared

testdir = "tests"


class ShutdownHandler(logging.Handler):
    """
    Trigger exit on errors.
    """

    def emit(self, record):
        logging.shutdown()
        sys.exit(1)


class DuplicateFilter(object):
    """
    Define a filter which keeps track of what was logged, and attach it to
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
# logging.getLogger().addFilter(DuplicateFilter())


__version__ = '0.9.0'
