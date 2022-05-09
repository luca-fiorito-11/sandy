# -*- coding: utf-8 -*-
"""
This module contains a single public function:

    * `get_df`

This function checks that a section defined by (mat, mf, mt) remains unchanged
after reading and writing it using the scripts inside the `sections` folder.

@author: abengoec
"""

import sandy
import pandas as pd
import io


def get_df_from_file(file, delimiter="?"):
    """
    Get the dataframe of a `endf6` file.

    Parameters
    ----------
    file : `str`
        Name of the `endf6` file.
    delimiter : `str`, optional
        Delimeter of the different lines. The default is "?".

    Returns
    -------
    `pd.DataFrame`
        Dataframe with the information of the file.

    """
    if isinstance(file, io.StringIO):
        text = file.read()
    else:
        with open(file) as f:
            text = f.read()

    def foo(x):
        return sandy.shared.add_delimiter_every_n_characters(
            x[:66],
            11,
            delimiter=delimiter,
        )
    newtext = "\n".join(map(foo, text.splitlines())).replace('"', '*')
    df = pd.read_csv(
        io.StringIO(sandy.shared.add_exp_in_endf6_text(newtext)),
        delimiter=delimiter,
        na_filter=True,
        names=["C1", "C2", "L1", "L2", "N1", "N2"],
    )
    df = df.apply(pd.to_numeric, errors='coerce')
    return df.fillna(0)


def get_df_from_text(text, delimiter="?"):
    """
    Get the dataframe of a `endf6` file in a text. 

    Parameters
    ----------
    text : `str`
        Text with the data in `endf6` format.
    delimiter : `str`, optional
        Delimeter of the different lines. The default is "?".

    Returns
    -------
    `pd.DataFrame`
        Dataframe with the information of the text introduced as input.

    Examples
    -------
    XS:
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'xs', 922350)
        >>> mt = 2
        >>> mat = 9228

        Obtain section information and replace blanks with zeros:
        >>> section = endf6._get_section_df(mat, 3, mt).apply(pd.to_numeric, errors='coerce').fillna(0)

        Get the same section using read and write functions:
        >>> read_xs = sandy.sections.mf3.read_mf3(endf6, mat, mt)
        >>> write_xs = sandy.sections.mf3.write_mf3(read_xs)
        >>> written_section = get_df_from_text(write_xs)

        Check the dataframes:
        >>> assert section.equals(written_section)
    """
    def foo(x):
        return sandy.shared.add_delimiter_every_n_characters(
            x[:66],
            11,
            delimiter=delimiter,
        )
    newtext = "\n".join(map(foo, text.splitlines())).replace('"', '*')
    df = pd.read_csv(
        io.StringIO(sandy.shared.add_exp_in_endf6_text(newtext)),
        delimiter=delimiter,
        na_filter=True,
        names=["C1", "C2", "L1", "L2", "N1", "N2"],
    )
    df = df.apply(pd.to_numeric, errors='coerce')
    return df.fillna(0)
