import logging
import sys

from .formats import read_formatted_file
from .settings import *
from .constants import *
from .decay import *
from .energy_grids import *
from .errorr import *
from .groupr import *
from .fy import *
from .h5file import *
from .libraries import *
from .pert import *
from .pfns import *
from .processing import *
from sandy.zam import *
from .njoy import *
from .sections import *
from .shared import *
from .utils import *
from .core import *
from .sampling import *
import sandy.mcnp
import sandy.aleph2
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
