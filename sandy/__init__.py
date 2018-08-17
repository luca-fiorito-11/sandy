import logging
import sys

from sandy.formats import *

class ShutdownHandler(logging.Handler):
    def emit(self, record):
        logging.shutdown()
        sys.exit(1)

FORMAT = '%(levelname)s:  %(message)s'
logging.basicConfig(format=FORMAT)
logging.getLogger().setLevel(logging.INFO)
logging.getLogger().addHandler(ShutdownHandler(level=40))