import pandas as pd
import pytest


collect_ignore = ["setup.py"]


@pytest.fixture(autouse=True, scope='session')
def pandas_terminal_width():
    pd.options.display.float_format = '{:.5e}'.format
    pd.set_option('display.max_columns', 1000)
    pd.set_option('display.width', 1000)
