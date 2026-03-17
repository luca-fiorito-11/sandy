"""
This module contains a single public function:

    * `read_mf35`

Function `read_mf35` reads a MF35/MT section from a string and produces a
content object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.
"""

__author__ = "Aitor Bengoechea"


def read_mf35(tape, mat, mt):
    """
    Parse MAT/MF=35/MT section from `sandy.Endf6` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number
    mt : `int`
        MT number

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------

    >>> import sandy
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 922380, local=True)
    >>> out = sandy.sections.mf35.read_mf35(tape, mat=9237, mt=18)
    >>> got = out["SUB"][1]["FKK"][0:15]
    >>> expected = [9.63626e-37, 6.16874e-36, 8.95968e-36, 3.04725e-35, 1.95073e-34, 2.8333e-34, 9.63626e-34,
    ...             6.16874e-33, 8.95968e-33, 3.04725e-32, 1.95073e-31, 2.8333e-31, 9.63625e-31, 6.16871e-30, 8.95957e-30]
    >>> assert got == expected
    """
    # --- IMPORT
    from ..records import read_cont_fast, read_list_fast

    mf = 35
    records = tape._get_section_records(mat, mf, mt)
    out = {"MAT": mat,
           "MF": mf,
           "MT": mt}
    i = 0
    C, i = read_cont_fast(records, i)
    out.update({"ZA": C.C1,
                "AWR": C.C2,
                "NK": C.N1, # Number of subsections
                "SUB": {}})
    for k in range(out["NK"]):
        L, i = read_list_fast(records, i)
        D = {"ELO": L.C1,  # Lowest incident neutron energy for this subsection
             "EHI": L.C2,  # Highest incident neutron energy for this subsection
             "LS": L.L1,  # Flago to indicate if the covariance matrix is symmetric
             "LB": L.L2,  # Flag to indicate if the covariance matrix is given in absolute or relative terms
             "NE": L.N2,  # Number of entries in the array containing outgoing particle energies
             "EK": L.B[:L.N2],  # Array containing outgoing particle energies
             "FKK": L.B[L.N2:]}  # Covariance matrix ordered by rows and starting from the diagonal term
        out["SUB"].update({k+1: D})
    return out
