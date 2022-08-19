import re
from io import StringIO

import pandas as pd


__author__ = "Luca Fiorito"

__all__ = [
        "get_keff",
        "get_table140",
        ]


def get_keff(file):
    with open(file, 'r') as f:
        text = f.read()
    PATTERN = "final estimate.*= (?P<keff>[\.0-9]+).*?(?P<stdev>[\.0-9]+)"
    match = re.search(PATTERN, text)
    keff = float(match.group("keff"))
    std = float(match.group("stdev"))
    return {"KEFF": keff, "ERR": std}


def get_table140(file):
    print(f"reading file '{file}'...")
    with open(file, 'r') as f:
        text = f.read()
    PATTERN = "(?:^1neutron.*table 140\n)(?P<table>(?:.*\n)+?)(?P<end>^\s+total\s{20})"
    match = re.search(PATTERN, text, re.MULTILINE).group("table")
    widths=[10, 9, 11, 9,] + [12] * 6
    dtypes = [int, int, str, float, int, float, float, float, float, int]
    names = pd.read_fwf(
        StringIO(match),
        widths=widths,
        skip_blank_lines=True,
        nrows=2,
        header=None,
        ).fillna("").apply(lambda col: " ".join(col).strip())
    df = pd.read_fwf(
        StringIO(match),
        widths=widths,
        skip_blank_lines=True,
        skiprows=3,
        names=names,
    ).ffill()
    return df.astype(dict(zip(names, dtypes)))