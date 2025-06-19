import logging
import sys

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


class ConditionalFormatter(logging.Formatter):
    """Change format dynamically based on log level."""
    
    FORMATS = {
        logging.INFO: logging.Formatter('%(message)s'),  # No level prefix for INFO
        'default': logging.Formatter('%(levelname)s:  %(message)s'),  # Default format
    }
    
    def format(self, record):
        formatter = self.FORMATS.get(record.levelno, self.FORMATS['default'])
        return formatter.format(record)

# Setup logging
logger = logging.getLogger()
handler = logging.StreamHandler()
handler.setFormatter(ConditionalFormatter())  # Apply the custom formatter
logger.addHandler(handler)
logger.setLevel(logging.INFO)


__version__ = '1.1.0'


# Must import submodules after __version__ and everything above.
from .constants import *
from .cov import *
from .decay import *
from .endf6 import *
from .energy_grids import *
from .errorr import *
from .fy import *
from .gendf import *
from .gls import *
from .libraries import *
from .lpc import *
from .pert import *
from .edistr import *
from .njoy import *
from .records import *
from .samples import *
from .sections import *
from .settings import *
from .shared import *
from .tools import *
from .tsl import *
from .utils import *
from .zam import *
from .sampling import *
from .spectra import *
from .xs import *

# These are folders
from . import mcnp
from . import aleph2
