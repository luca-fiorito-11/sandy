import pandas as pd
import pytest
import logging
import os

@pytest.fixture(autouse=True, scope='session')
def configure_environment():
    # Pandas display settings
    pd.options.display.float_format = '{:.5e}'.format
    pd.set_option('display.max_columns', 1000)
    pd.set_option('display.width', 1000)

    # Determine log level from environment variable or default to INFO
    log_level = os.getenv("PYTEST_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, log_level, logging.INFO)

    # Override root logger configuration
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Remove existing handlers to avoid duplicates
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Create a new handler with conditional formatting
    handler = logging.StreamHandler()
    is_ci = os.getenv("CI", "false").lower() == "true"
    if is_ci:
        formatter = logging.Formatter("%(levelname)s: %(message)s")
    else:
        formatter = logging.Formatter("[%(asctime)s] %(levelname)s in %(module)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)
