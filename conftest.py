import pandas as pd
import pytest


collect_ignore = ["setup.py"]

pd.options.display.float_format = '{:.5e}'.format


@pytest.fixture(autouse=True, scope='session')
def pandas_terminal_width():
    pd.set_option('display.max_columns', 1000)
