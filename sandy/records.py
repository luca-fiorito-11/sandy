from typing import Tuple, List
import math

__author__ = "Luca Fiorito"


line_pattern = "{:<66}{:4d}{:2d}{:3d}{:5d}"


"""
Helper functions used in ``read_cont``.
"""
def _is_empty(val):
    """Return True for None, NaN, or whitespace-only strings."""
    if val is None:
        return True
    if isinstance(val, float) and math.isnan(val):
        return True
    if isinstance(val, str) and val.strip() == "":
        return True
    return False


def _to_float(val):
    """Convert to float, treating empty values as 0.0."""
    if _is_empty(val):
        return 0.0

    s = str(val).strip()
    try:
        return float(s)  # allow exponent, trailing decimal, etc.
    except Exception:
        raise


def _to_int(val):
    """
    Convert to int, treating empty values as 0 and normalizing '3.0' -> 3.
    
    >>> import sandy, pytest
    >>> assert sandy.records._to_int(3.) == 3
    >>> assert sandy.records._to_int(3e0) == 3
    >>> assert sandy.records._to_int(3.0) == 3
    >>> assert sandy.records._to_int("3.0") == 3
    >>> assert sandy.records._to_int("       3.0") == 3
    >>> with pytest.raises(Exception):
    ...    sandy.records._to_int(3.2)
    """
    if _is_empty(val):
        return 0

    f = _to_float(val)

    # Check: is it mathematically an integer?
    if not f.is_integer():
        raise

    return int(f)



def write_cont(C1, C2, L1, L2, N1, N2):
    """
    Write ENDF-6 **cont** record.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string

    Warns
    -----
    This function will produce strings longer than 66 characters if integers
    `> 99999999999` are given.
    """
    integers = f"{L1:11d}{L2:11d}{N1:11d}{N2:11d}"
    line = write_float(C1) + write_float(C2) + integers
    return [line]


def write_text(text):
    """
    Write ENDF-6 `TEXT` record in formatted fortran.

    Returns
    -------
    `str`
        list of 66-characters-long ENDF-6 formatted string

    """
    line = f"{text[:66]:66}"
    return [line]


def write_integer_list(lst):
    """
    Write list of integers into ENDF-6 format.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string
    """
    from .utils import grouper

    itr = grouper(map("{:11d}".format, lst), 6, fillvalue=" "*11)

    return ["".join(vals) for vals in itr]


def write_float_list(lst):
    """
    Write list of floats into ENDF-6 format.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string
    """
    from .utils import grouper

    itr = grouper(map(write_float, lst), 6, fillvalue=" "*11)

    return ["".join(vals) for vals in itr]



def write_tab2(C1, C2, L1, L2, N2, NBT, INT):
    """
    Write ENDF-6 **tab2** record.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string
    """
    from .utils import interwine_lists

    N1 = len(NBT)
    lines = write_cont(C1, C2, L1, L2, N1, N2)
    lines += write_integer_list(interwine_lists(NBT, INT))
    return lines



def write_tab1(C1, C2, L1, L2, NBT, INT, x, y):
    """
    Write ENDF-6 **tab1** record.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string
    """
    from .utils import interwine_lists

    N2 = len(x)
    lines = write_tab2(C1, C2, L1, L2, N2, NBT, INT)
    lines += write_float_list(interwine_lists(x, y))

    return lines



def get_records_from_text(
        text: str,
        ) -> List[Tuple[str, str, str, str, str, str]]:
    """
    Parse an ENDF-6 text block into a list of fixed-width ENDF-6 records.

    Each ENDF-6 data line contains 66 characters of numerical content,
    organized into six 11-character fields (C1, C2, L1, L2, N1, N2).
    Any trailing MAT/MF/MT identifiers (in columns 67–80) are ignored.

    This function splits the text into lines and extracts the six
    fixed-width fields from each line using ``split_endf6_line``.

    Parameters
    ----------
    text : str
        Multiline string containing ENDF-6 formatted records.

    Returns
    -------
    list of tuple(str, str, str, str, str, str)
        A list where each element is a 6-tuple:
        ``(C1, C2, L1, L2, N1, N2)``,
        all unprocessed fixed-width string slices of 11 characters each.

    Raises
    ------
    TypeError
        If ``text`` is not a string.
    ValueError
        If at least one line is shorter than 66 characters
        (caught inside ``split_endf6_line``).

    Examples
    --------
    
    First 5 lines of text from the JEFF-3.3 decay data file of U-235.

    >>> text = '''9.223500+4 2.330250+2          0          0          0          63542 8457    1
    ...  2.22102+16 1.57788+13          0          0          6          03542 8457    2
    ...  5.067170+4 4.291630+3 1.636160+5 1.708010+3 4.464600+6 1.632550+53542 8457    3
    ...  3.500000+0-1.000000+0          0          0         12          23542 8457    4
    ...  4.000000+0 0.000000+0 4.678700+6 7.000000+2 1.000000+0 1.000000-43542 8457    5'''
    >>> expected = [
    ...    ('9.223500+4 ', '2.330250+2 ', '         0 ', '         0 ', '         0 ', '         63'),
    ...    (' 2.22102+16', ' 1.57788+13', '          0', '          0', '          6', '          0'),
    ...    (' 5.067170+4', ' 4.291630+3', ' 1.636160+5', ' 1.708010+3', ' 4.464600+6', ' 1.632550+5'),
    ...    (' 3.500000+0', '-1.000000+0', '          0', '          0', '         12', '          2'),
    ...    (' 4.000000+0', ' 0.000000+0', ' 4.678700+6', ' 7.000000+2', ' 1.000000+0', ' 1.000000-4')]
    >>> got = get_records_from_text(text)
    >>> assert got == expected
    """
    if not isinstance(text, str):
        raise TypeError(f"get_records_from_text expected a string, got {type(text).__name__}")

    records = []

    for line in text.splitlines():
        C1, C2, L1, L2, N1, N2 = split_endf6_line(line)
        records.append((C1, C2, L1, L2, N1, N2))

    return records


def split_endf6_line(
        line: str,
        ) -> Tuple[str, str, str, str, str, str]:
    """
    Split an ENDF-6 data line into its six fixed-width fields (C1, C2, L1, L2, N1, N2).

    Parameters
    ----------
    line : str
        A full ENDF-6 line. The function extracts the first 66 characters,
        which correspond to six fixed-width 11-character fields.
        Any trailing MAT/MF/MT identifiers (columns 67–80) are ignored.

    Returns
    -------
    tuple of str
        A 6-tuple of raw string fields (C1, C2, L1, L2, N1, N2),
        each exactly 11 characters long as defined by the ENDF-6 format.

    Raises
    ------
    TypeError
        If `line` is not a string.
    ValueError
        If the line is shorter than 66 characters.

    Notes
    -----
    ENDF-6 lines have a 66-character data region (columns 1–66),
    followed by optional MAT/MF/MT numbering (columns 67–80).
    This function discards everything beyond index 66.

    Examples
    --------
    Standard 66-charachter line.

    >>> line = ' 9.223500+4 2.330250+2          0          0          0          6'
    >>> expected = (' 9.223500+4', ' 2.330250+2', '          0', '          0', '          0', '          6')
    >>> got = split_endf6_line(line)
    >>> assert got == expected

    If MAT/MF/MT are given in the line, they are not reported in outupt.
    
    >>> line = ' 9.223500+4 2.330250+2          0          0          0          63542 8457    1'
    >>> expected = (' 9.223500+4', ' 2.330250+2', '          0', '          0', '          0', '          6')
    >>> got = split_endf6_line(line)
    >>> assert got == expected
    
    Fail if non-string is given.

    >>> import pytest
    >>> with pytest.raises(Exception):
    ...    split_endf6_line(3)

    Fail if string does not have 66 chaacters.

    >>> import pytest
    >>> with pytest.raises(Exception):
    ...    split_endf6_line("dsf d")
    """
    # ---- Type check ---
    if not isinstance(line, str):
        raise TypeError(f"split_line expected a string, got {type(line).__name__}")

    # ---- Length check ---
    if len(line) < 66:
        raise ValueError(
            f"ENDF-6 line must be at least 66 characters, got {len(line)}: {line!r}"
        )

    # ---- Extract the fixed-width 66-character data region ---
    data = line[:66]

    # --- Slice into six 11-character fields ---
    C1 = data[0:11]
    C2 = data[11:22]
    L1 = data[22:33]
    L2 = data[33:44]
    N1 = data[44:55]
    N2 = data[55:66]

    return C1, C2, L1, L2, N1, N2



class TextRecord:
    """
    Lightweight container for an ENDF-6 TEXT record.

    A TEXT record consists of 66 characters taken from the six
    fixed-width fields (C1, C2, L1, L2, N1, N2) of a single ENDF-6 line.
    """
    __slots__ = ("HL",)
    
def read_text_fast(
        records: List[Tuple[str, str, str, str, str, str]],
        i: int,
    ) -> Tuple[TextRecord, int]:
    """
    Read a single ENDF-6 TEXT record from fixed-width string fields.

    ENDF-6 TEXT records contain descriptive formatted text encoded
    in the six 11-character fields (66 characters total). This function
    reconstructs that 66-character string exactly as stored.

    Parameters
    ----------
    records : list of tuple(str, str, str, str, str, str)
        Fixed-width ENDF-6 fields extracted using `split_endf6_line`.

    i : int
        Index of the TEXT record to read.

    Returns
    -------
    (TextRecord, int)
        A tuple containing:
        - the parsed TEXT record
        - the updated index (i + 1)

    Raises
    ------
    TypeError
        If `records` is not a list or `i` is not an int.
    IndexError
        If `i` is out of range.
    ValueError
        If the record does not contain exactly six fields.


    Examples
    --------
    >>> text = get_records_from_text(
    ...     'DECAY DATA FOR URANIUM-235                             9228 8457    1'
    ... )
    >>> rec, j = read_text_fast(text, 0)
    >>> rec.HL.startswith("DECAY DATA")
    True
    >>> j
    1
    """

    if not isinstance(records, list):
        raise TypeError(f"records must be a list, got {type(records).__name__}")
    if not isinstance(i, int):
        raise TypeError(f"i must be int, got {type(i).__name__}")

    try:
        fields = records[i]
    except Exception:
        raise IndexError(f"Index {i} out of range for ENDF records")

    if not isinstance(fields, tuple) or len(fields) != 6:
        raise ValueError(
            f"TEXT record must contain 6 fixed-width fields, got: {fields}"
        )

    # Construct the 66-character TEXT field
    C1, C2, L1, L2, N1, N2 = fields
    HL = C1 + C2 + L1 + L2 + N1 + N2

    rec = TextRecord()
    rec.HL = HL

    return rec, i + 1



class ContRecord:
    """
    Lightweight container for an ENDF-6 CONT record.
    Fields correspond to the six fixed-width floats/ints:
    C1, C2, L1, L2, N1, N2.
    """
    __slots__ = (
        "C1","C2","L1","L2","N1","N2"
        )

def read_cont_fast(
        records: List[Tuple[str, str, str, str, str, str]],
        i: int,
        ) -> Tuple[ContRecord, int]:
    """
    Parse a single ENDF-6 CONT record from a list of fixed-width string fields.

    Parameters
    ----------
    records : list of tuple(str, str, str, str, str, str)
        Each element of the list is a 6-tuple of 11-character strings
        (C1, C2, L1, L2, N1, N2), obtained from splitting ENDF-6 lines
        via `split_endf6_line` or similar.

    i : int
        Index of the CONT record to read. This function consumes exactly
        one ENDF-6 line and returns the updated pointer.

    Returns
    -------
    (ContRecord, int)
        A tuple containing:
        - the parsed `ContRecord` instance
        - the updated index `i + 1`

    Raises
    ------
    TypeError
        If `records` is not a list, or fields are malformed.
    IndexError
        If `i` is out of range.
    ValueError
        If the record does not contain exactly six fields.

    Notes
    -----
    A CONT record contains:
        C1 : float
        C2 : float
        L1 : int
        L2 : int
        N1 : int
        N2 : int

    The caller is responsible for ensuring that the fields are
    unprocessed fixed-width ENDF-6 strings and will be converted here
    using `parse_endf_float` and `parse_endf_int`.

    - Empty strings (including whitespace-only) and NaN are treated as 0 for all fields.
    - If something is wrong, an error is reported with the line number (1-based) and the line is printed.

    Examples
    --------
    Test some features.

    >>> import sandy, pytest

    Minimal test on the JEFF-3.3 decay data for U-235.

    >>> records = [
    ...    (' 9.223500+4', ' 2.330250+2', '          0', '          0', '          0', '          6')
    ... ]
    >>> rec, new_i = sandy.records.read_cont_fast(records, 0)
    >>> assert (rec.C1, rec.C2, rec.L1, rec.L2, rec.N1, rec.N2) == (92235.0, 233.025, 0, 0, 0, 6)
    >>> assert new_i == 1

    This happens often at the end of the intro section MF=1 MT=451.
    Empty sections are converted to zeros.

    >>> line = "                                1        451         29          01960 1451   28"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(1960, 1, 451)
    >>> out, ipos = sandy.records.read_cont_fast(records, 0)
    >>> assert out.C1 == 0 and out.C2 == 0 and out.L1 == 1 and out.L2 == 451 and out.N1 == 29 and out.N2 == 0
    >>> assert ipos == 1

    This is a line found in an intro section MF=1 MT=451, but it is also the standard way.
    
    >>> line = " 0.00000+00 0.00000+00          8        457         41          02007 1451   51"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(2007, 1, 451)
    >>> out, ipos = sandy.records.read_cont_fast(records, 0)
    >>> assert out.C1 == 0 and out.C2 == 0 and out.L1 == 8 and out.L2 == 457 and out.N1 == 41 and out.N2 == 0
    
    An error 'xxx' is inserted.

    >>> line = "  xxx                           1        451         29          01960 1451   28"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(1960, 1, 451)
    >>> with pytest.raises(Exception):
    ...    sandy.records.read_cont_fast(records, 0)
    
    A number '2' is shifted but still recognized.

    >>> line = "  2                             1        451         29          01960 1451   28"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(1960, 1, 451)
    >>> out, ipos = sandy.records.read_cont_fast(records, 0)
    >>> assert out.C1 == 2.0 and out.C2 == 0 and out.L1 == 1 and out.L2 == 451 and out.N1 == 29 and out.N2 == 0
    
    Let's have floats everywhere, but N1, N2, L1 and L2 can be coverted to ``int``.

    >>> line = " 1.00000+00 1.00000+00 1.00000+00 1.00000+00 1.00000+00 1.00000+002007 1451   51"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(2007, 1, 451)
    >>> out, ipos = sandy.records.read_cont_fast(records, 0)
    >>> assert out.C1 == out.C2 == out.L1 == out.L2 == out.N1 == out.N2 == 1

    Let's have floats everywhere, but N1, N2, L1 and L2 cannot be coverted to ``int``.
    An error should pop up.

    >>> line = " 1.00000+00 1.00000+00 1.00000+00 1.00000+00 1.00000+00 1.10000+002007 1451   51"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(2007, 1, 451)
    >>> with pytest.raises(Exception):
    ...    sandy.records.read_cont_fast(records, 0)

    >>> line = "                                                                  2007 1451   51"
    >>> records = sandy.Endf6.from_text(line)._get_section_records(2007, 1, 451)
    >>> out, ipos = sandy.records.read_cont_fast(records, 0)
    >>> assert out.C1 == out.C2 == out.L1 == out.L2 == out.N1 == out.N2 == 0
    """

    C1s, C2s, L1s, L2s, N1s, N2s = records[i]
    rec = ContRecord()
    rec.C1 = parse_endf_float(C1s)
    rec.C2 = parse_endf_float(C2s)
    rec.L1 = parse_endf_int(L1s)
    rec.L2 = parse_endf_int(L2s)
    rec.N1 = parse_endf_int(N1s)
    rec.N2 = parse_endf_int(N2s)
    return rec, i+1



class ListRecord:
    """
    Lightweight container for an ENDF-6 LIST record.
    A LIST record consists of:
        - a CONT header: C1, C2, L1, L2, NPL, N2
        - a list of N1 float values stored in B
    """
    __slots__ = (
        "C1", "C2", "L1", "L2", "NPL", "N2", "B"
        )

def read_list_fast(
        records: List[Tuple[str, str, str, str, str, str]],
        i: int,
        ) -> Tuple[ListRecord, int]:
    """
    Reads:
        - one CONT record using `read_cont_fast`
        - then N1 floating-point values, where N1 = C.N1

    Parameters
    ----------
    records : list of tuple(str, str, str, str, str, str)
        List of fixed-width ENDF-6 line fields produced by `split_endf6_line`.

    i : int
        Index of the CONT record to read.

    Returns
    -------
    (ListRecord, int)
        The parsed LIST record, and the updated index.

    Raises
    ------
    TypeError
        If input types are incorrect.
    IndexError
        If `i` is out of range.
    ValueError
        If the fields are malformed.

    Notes
    -----
    This function implements the generic ENDF LIST rule:
        LIST length = C.N1

    Special formats (TAB1, TAB2, MF=8 energy LIST, etc.)
    that require LIST lengths other than C.N1 must use
    `read_list_with_size_fast()` instead.

    Examples
    --------
    Minimal test, reading half-live part of U-235 decay file from JEFF-3.3.

    >>> text = [
    ...    ' 2.22102+16 1.57788+13          0          0          6          0',
    ...    ' 5.067170+4 4.291630+3 1.636160+5 1.708010+3 4.464600+6 1.632550+5',
    ... ]
    >>> records = get_records_from_text('\\n'.join(text))
    >>> i = 0
    >>> L, j = read_list_fast(records, i)
    >>> expected = (2.22102e+16, 1.57788e+13, 0, 0, 6, 0)
    >>> got = (L.C1, L.C2, L.L1, L.L2, L.NPL, L.N2)
    >>> assert expected == got
    >>> assert j == 2
    >>> expected = [5.067170e+4, 4.291630e+3, 1.636160e+5, 1.708010e+3, 4.464600e+6, 1.632550e+5]
    >>> got = L.B
    >>> assert expected == got
    """

    # --- Type checks ---
    if not isinstance(records, list):
        raise TypeError(f"`records` must be a list, got {type(records).__name__}")
    if not isinstance(i, int):
        raise TypeError(f"`i` must be int, got {type(i).__name__}")

    # --- 1) Read CONT record ---
    C, i = read_cont_fast(records, i)

    # --- 2) LIST size = C.N1 ---
    size = C.N1

    # --- 3) Read LIST values ---
    B, i = read_list_with_size_fast(records, i, size)

    # --- 4) Build record ---
    rec = ListRecord()
    rec.C1 = C.C1
    rec.C2 = C.C2
    rec.L1 = C.L1
    rec.L2 = C.L2
    rec.NPL = C.N1
    rec.N2 = C.N2
    rec.B = B

    return rec, i

def read_list_with_size_fast(
    records: List[Tuple[str, str, str, str, str, str]],
    i: int,
    size: int,
) -> Tuple[List[float], int]:
    """
    Read a LIST of `size` floating-point values from ENDF-6 records.

    This function does NOT read a CONT record.

    Parameters
    ----------
    records : list of tuple(str, ..., str)
        List of fixed-width ENDF field tuples.

    i : int
        Index of the first value row.

    size : int
        Number of float values to read.

    Returns
    -------
    (list of float, int)
        The parsed float values, and the updated pointer.

    Raises
    ------
    TypeError
        If inputs have wrong types.
    IndexError
        If insufficient rows exist.
    ValueError
        If size is negative.

    Notes
    -----
    LIST values are stored 6 per ENDF-6 line.
    ENDF-6 standard: number of rows = ceil(size / 6).
    """

    from itertools import chain
    import math

    # --- Input checks ---
    if not isinstance(records, list):
        raise TypeError(f"`records` must be a list, got {type(records).__name__}")
    if not isinstance(i, int):
        raise TypeError(f"`i` must be int, got {type(i).__name__}")
    if not isinstance(size, int):
        raise TypeError(f"`size` must be int, got {type(size).__name__}")
    if size < 0:
        raise ValueError(f"`size` must be >= 0, got {size}")

    # --- Determine how many rows hold value fields ---
    iadd = math.ceil(size / 6)

    # --- Extract rows ---
    rows = records[i : i + iadd]
    if len(rows) < iadd:
        raise IndexError(
            f"Not enough ENDF lines to read LIST of size {size}: "
            f"needed {iadd}, only {len(rows)} available."
        )

    # --- Flatten rows and convert ---
    flat = list(chain.from_iterable(rows))
    B = [parse_endf_float(v) for v in flat[:size]]

    # --- Advance pointer ---
    i += iadd

    return B, i



class Tab2Record:
    """
    Lightweight container for an ENDF-6 TAB2 record.

    A TAB2 record controls a 1D or 2D tabulation structure.
    It stores:

    C1, C2 : float
        Context-dependent parameters (e.g., ZA, AWR, Q-values).
    L1, L2 : int
        Format-dependent modifiers or flags.
    NR : int
        Number of interpolation ranges.
    NZ : int
        Number of subfunctions that follow (e.g. number of TAB1/LIST blocks).
        In TAB1 usage, NZ is the number of (x,y) pairs.
        In TAB2 usage, NZ is the number of z-values or secondary tabs.

    NBT : list of int
        Upper boundaries of interpolation intervals (length = NR).
    INT : list of int
        Interpolation law identifiers for each interval (length = NR).

    Notes
    -----
    TAB2 is always followed by a LIST of length 2*NR, containing:
        (NBT_1, INT_1), (NBT_2, INT_2), ..., (NBT_NR, INT_NR)

    This object is a passive container. All parsing logic is handled
    by `read_tab2_fast`.
    """
    __slots__ = (
        "C1", "C2", "L1", "L2", "NR", "NZ", "NBT", "INT"
        )


def read_tab2_fast(
        records: List[Tuple[str, str, str, str, str, str]],
        i: int,
        ) -> Tuple[Tab2Record, int]:
    """
    Reads:
        1. A CONT record: (C1, C2, L1, L2, NR, NZ)
        2. A LIST record of size 2*NR containing pairs of (NBT, INT)

    Parameters
    ----------
    records : list of tuple(str, str, str, str, str, str)
        List of fixed-width ENDF-6 fields obtained from `split_endf6_line`.

    i : int
        Index of the first line of the TAB2 record.

    Returns
    -------
    (Tab2Record, int)
        The parsed TAB2 record and the updated index.

    Raises
    ------
    TypeError
        If records or i have invalid types.
    IndexError
        If there are insufficient records to read TAB2.
    ValueError
        If malformed data are detected.

    Notes
    -----
    This function is equivalent to:

        C, ipos = read_cont(...)
        L, ipos = _read_list(df, ipos, 2*C.N1)

    but implemented without DataFrames, making it 50–100x faster.

    Examples
    --------
    Typical usage on a section taken from the U235 decay file for JEFF-3.3.
    Same section used to test :func:`sandy.records.read_tab1_fast`.

    >>> text = [' 6.000000+0 0.000000+0          0          0          1         15',
    ...  '         15          1                                            ']
    >>> records = get_records_from_text("\\n".join(text))
    >>> i = 0
    >>> T, j = read_tab2_fast(records, i)
    >>> expected = 6.000000e+0, 0.000000e+0,          0,          0,          1,         15
    >>> got = T.C1, T.C2, T.L1, T.L2, T.NR, T.NZ
    >>> assert expected == got
    >>> assert j == 2
    >>> assert T.NBT == [15]
    >>> assert T.INT == [1]
    """
    # ---- Input checks ----
    if not isinstance(records, list):
        raise TypeError(f"`records` must be list, got {type(records).__name__}")
    if not isinstance(i, int):
        raise TypeError(f"`i` must be int, got {type(i).__name__}")

    # ---- 1. Read CONT ----
    C, i = read_cont_fast(records, i)

    NR = C.N1
    NZ = C.N2

    # ---- 2. Read LIST of length 2*NR ----
    size = 2 * NR
    B, i = read_list_with_size_fast(records, i, size)

    # B already contains floats → convert to integers
    NBT = [int(v) for v in B[0::2]]
    INT = [int(v) for v in B[1::2]]

    # ---- 3. Build TAB2 record ----
    rec = Tab2Record()
    rec.C1 = C.C1
    rec.C2 = C.C2
    rec.L1 = C.L1
    rec.L2 = C.L2
    rec.NR = NR
    rec.NZ = NZ
    rec.NBT = NBT
    rec.INT = INT

    return rec, i



class Tab1Record:
    """
    Lightweight container for an ENDF-6 TAB1 record.

    A TAB1 record defines a one-dimensional tabulated function y(x)
    together with interpolation rules. It consists of:

    - A controlling TAB2 record containing:
        C1, C2 : float
            (context-dependent parameters)
        L1, L2 : int
            (context-dependent flags)
        NR     : int
            Number of interpolation ranges
        NP     : int
            Number of (x, y) tabulated pairs

    - Interpolation tables:
        NBT : list of int
            Upper bounds of interpolation intervals
        INT : list of int
            Interpolation law identifiers for each interval

    - Tabulated function values:
        x : list of float
            The independent variable values x(n)
        y : list of float
            The corresponding dependent variable values y(n)

    Notes
    -----
    TAB1 is used throughout many MF sections. This object only stores the parsed content;
    all pointer logic is handled by `read_tab1_fast`.
    """
    __slots__ = (
        "C1","C2","L1","L2","NR","NP","NBT","INT","x","y"
        )

def read_tab1_fast(
        records: List[Tuple[str, ...]], 
        i: int
        ) -> Tuple[Tab1Record, int]:
    """
    Fast, fully compatible replacement for `sandy.records.read_tab1`.

    Reads an ENDF-6 TAB1 record, which consists of:

        1. A TAB2 record:
            - one CONT line
            - one LIST of length 2*NR (pairs of NBT and INT)

        2. A LIST of length 2*NP (NP = NZ from TAB2)
           containing x0, y0, x1, y1, ..., x(NP-1), y(NP-1).

    Parameters
    ----------
    records : list of tuple(str, str, str, str, str, str)
        Fixed-width ENDF-6 field tuples extracted by `split_endf6_line`.

    i : int
        Index of the first line of the TAB1 record.

    Returns
    -------
    (Tab1Record, int)
        The parsed TAB1 record and the updated index.

    Raises
    ------
    TypeError
        If inputs are of incorrect type.
    IndexError
        If insufficient records exist to complete the TAB1 structure.
    ValueError
        If ENDF field content is malformed.

    Examples
    --------
    Typical usage on a section taken from the U235 decay file for JEFF-3.3.

    >>> text = [' 6.000000+0 0.000000+0          0          0          1         15',
    ...  '         15          1                                            ',
    ...  ' 0.000000+0 2.857140-7 1.400000+5 9.246810-7 3.000000+5 1.021120-6',
    ...  ' 5.000000+5 8.828020-7 7.000000+5 4.973310-7 1.000000+6 2.594850-7',
    ...  ' 1.500000+6 1.201260-7 2.000000+6 7.232010-8 2.500000+6 4.559450-8',
    ...  ' 3.000000+6 2.104910-8 4.000000+6 7.927030-9 5.000000+6 2.955160-9',
    ...  ' 6.000000+6 1.162160-9 7.000000+6 7.05868-11 1.000000+7 0.000000+0']
    >>> records = get_records_from_text("\\n".join(text))
    >>> i = 0
    >>> T, j = read_tab1_fast(records, i)
    >>> expected = 6.000000e+0, 0.000000e+0,          0,          0,          1,         15
    >>> got = T.C1, T.C2, T.L1, T.L2, T.NR, T.NP
    >>> assert expected == got
    >>> assert j == 7
    >>> assert T.NBT == [15]
    >>> assert T.INT == [1]
    >>> expected = [0.000000e+0, 1.400000e+5, 3.000000e+5, 5.000000e+5, 7.000000e+5, 1.000000e+6, 1.500000e+6,
    ...             2.000000e+6, 2.500000e+6, 3.000000e+6, 4.000000e+6, 5.000000e+6, 6.000000e+6, 7.000000e+6, 1.000000e+7]
    >>> assert T.x == expected
    >>> expected = [2.857140e-7, 9.246810e-7, 1.021120e-6, 8.828020e-7, 4.973310e-7, 2.594850e-7, 1.201260e-7,
    ...             7.232010e-8, 4.559450e-8, 2.104910e-8, 7.927030e-9, 2.955160e-9, 1.162160e-9, 7.05868e-11, 0.000000e+0]
    >>> assert T.y == expected
    """
    # --- Lightweight input checks ---
    if not isinstance(records, list):
        raise TypeError(f"`records` must be list, got {type(records).__name__}")
    if not isinstance(i, int):
        raise TypeError(f"`i` must be int, got {type(i).__name__}")

    # --- 1) Read TAB2 (CONT + LIST(2 NR)) ---
    T2, i = read_tab2_fast(records, i)

    NR = T2.NR     # number of interpolation ranges
    NP = T2.NZ     # number of (x,y) pairs

    # --- 2) Read LIST of length 2*NP ---
    size = 2 * NP
    B, i = read_list_with_size_fast(records, i, size)

    # B already contains floats; just split them
    x = B[0::2]
    y = B[1::2]

    # --- 3) Build final record ---
    rec = Tab1Record()
    rec.C1 = T2.C1
    rec.C2 = T2.C2
    rec.L1 = T2.L1
    rec.L2 = T2.L2
    rec.NR = NR
    rec.NP = NP
    rec.NBT = T2.NBT
    rec.INT = T2.INT
    rec.x = x
    rec.y = y

    return rec, i



def parse_endf_float(x: str) -> float:
    """
    Convert ENDF-6 float fields into Python float.

    Handles two formats:
    - ENDF Fortran-style without 'E'  (e.g. ' 1.0000+1')
    - Standard scientific notation with 'E' (e.g. ' 1.0000E+01')

    Parameters
    ----------
    x : str
        ENDF-6 float field (11 chars), may use Fortran-style exponent
        without 'E', or standard 'E' notation.

    Returns
    -------
    float
        Parsed value. Blank fields return 0.0.

    Examples
    --------
    >>> parse_endf_float(" 1.000000+1")
    10.0
    >>> parse_endf_float("-3.21540-3 ")
    -0.0032154
    >>> parse_endf_float(" 1.0000E+00")
    1.0
    >>> parse_endf_float("           ")
    0.0
    >>> assert parse_endf_float(" 1.000000+1") == 10.0
    >>> assert abs(parse_endf_float("-3.21540-3 ") + 0.0032154) < 1e-12
    >>> assert parse_endf_float(" 1.0000E+00") == 1.0
    >>> assert parse_endf_float(" 5.000000E-01") == 0.5
    >>> assert parse_endf_float(" 1.0000e+00") == 1.0
    >>> assert parse_endf_float(" 5.000000e-01") == 0.5
    >>> assert parse_endf_float("           ") == 0.0
    >>> assert parse_endf_float(" 12.34") == 12.34
    >>> import pytest
    >>> with pytest.raises(Exception):
    ...    assert parse_endf_float("ABCDE")
    """
    if not isinstance(x, str):
        raise TypeError(f"parse_endf_float expected str, got {type(x).__name__}")

    s = x.strip()
    if not s:
        return 0.0
    
    # If number already uses E-format, use Python float directly.
    if 'E' in s or 'e' in s:
        return float(s)

    # ENDF-style format: find exponent sign after the first character.
    for i in range(1, len(s)):
        c = s[i]
        if c in "+-":
            # Insert "E" before the exponent
            return float(s[:i] + "E" + s[i:])

    # No exponent → normal float literal
    return float(s)



def parse_endf_int(x: str) -> int:
    """
    Convert an ENDF-6 integer field into Python int.

    ENDF-6 integer fields are stored as right-aligned numbers within
    an 11-character field. Blank or whitespace-only fields represent 0.

    Parameters
    ----------
    x : str
        A fixed-width ENDF-6 integer field.

    Returns
    -------
    int
        Parsed integer value.

    Raises
    ------
    TypeError
        If `x` is not a string.

    Notes
    -----
    This function is used in:
        - `read_cont_fast` (L1, L2, N1, N2)

    Examples
    --------
    >>> parse_endf_int('         12')
    12
    >>> parse_endf_int('         3.0')
    3
    >>> parse_endf_int('          ')
    0
    >>> import pytest
    >>> with pytest.raises(Exception):
    ...    parse_endf_int('         0.2')
    """
    x = parse_endf_float(x)
    if isinstance(x, float) and not x.is_integer():
        raise ValueError(f"{x} is not an integer value")

    return int(x)



def write_list(C1, C2, L1, L2, N2, B):
    """
    Write ENDF-6 **list** record.

    Outputs:
        - list of string
    """
    NPL = len(B)
    lines = write_cont(C1, C2, L1, L2, NPL, N2)
    lines += write_float_list(B)
    return lines


def write_float(x):
    """
    Converts a floating-point number to a formatted string representation.

    Parameters
    ----------
        x (float): The floating-point number to be converted.

    Returns
    -------
        str: The formatted string representation of the floating-point number.

    Examples
    --------

    >>> import sandy
    >>> assert sandy.records.write_float(2) == ' 2.00000000'
    >>> assert sandy.records.write_float(2e1) == ' 20.0000000'
    >>> assert sandy.records.write_float(2e2) == ' 200.000000'
    >>> assert sandy.records.write_float(2e3) == ' 2000.00000'
    >>> assert sandy.records.write_float(2e4) == ' 20000.0000'
    >>> assert sandy.records.write_float(2e5) == ' 200000.000'
    >>> assert sandy.records.write_float(2e6) == ' 2000000.00'
    >>> assert sandy.records.write_float(2e7) == ' 20000000.0'
    >>> assert sandy.records.write_float(2e8) == '  200000000'
    >>> assert sandy.records.write_float(0) == ' 0.00000000'
    >>> assert sandy.records.write_float(2e-3) == ' 2.000000-3'
    >>> assert sandy.records.write_float(2e-10) == ' 2.00000-10'
    >>> assert sandy.records.write_float(1-1e-8) == ' 1.000000+0'
    """
    if abs(x) >= 1e0 and abs(x) < 1e1:
        y = f"{x:11.8f}"
    elif abs(x) >= 1E1 and abs(x) < 1E2:
        y = f"{x:11.7f}"
    elif abs(x) >= 1E2 and abs(x) < 1E3:
        y = f"{x:11.6f}"
    elif abs(x) >= 1E3 and abs(x) < 1E4:
        y = f"{x:11.5f}"
    elif abs(x) >= 1E4 and abs(x) < 1E5:
        y = f"{x:11.4f}"
    elif abs(x) >= 1E5 and abs(x) < 1E6:
        y = f"{x:11.3f}"
    elif abs(x) >= 1E6 and abs(x) < 1E7:
        y = f"{x:11.2f}"
    elif abs(x) >= 1E7 and abs(x) < 1E8:
        y = f"{x:11.1f}"
    elif abs(x) >= 1E8 and abs(x) < 1E10:
        y = f"{x:11.0f}"
    elif x == 0:
        y = f"{x:11.8f}"
    elif abs(x) < 1E0 and abs(x) >= 1e-9:
        y = f"{x:13.6e}".replace("e-0", "-").replace("e+0", "+")
    else:
        y = f"{x:12.5e}".replace("e", "")
    return y


def write_int(x):
    """
    Examples
    --------
    
    >>> import sandy
    >>> sandy.records.write_int(10)
    '         10'

    >>> sandy.records.write_int(-1e5)
    '    -100000'

    >>> import pytest
    >>> with pytest.raises(ValueError): sandy.records.write_int(-1e10)
    """
    y = f"{int(x):>11d}"
    if len(y) > 11:
        raise ValueError(f"Integer '{y}' exceeds 11 characters")
    return y


def line_numbers(length):
    """
    Line number creator

    Parameters
    ----------
    length : 'int'
        Number of lines.

    Returns
    -------
    ilines : `list`
        List containing the number of each line.

    Examples
    --------
    Some tests.

    >>> assert max(line_numbers(1.0e6)) == 99999
    >>> assert min(line_numbers(1.0e6+1)) == 1
    >>> assert max(line_numbers(1.0e4+1)) == 10001
    >>> assert len(line_numbers(1.0e6)) == 1000000
    >>> assert len(line_numbers(1.0e6+1)) == 1000001
    """
    import numpy as np

    iend = 1 + length
    ilines = np.tile(np.arange(1, 1e5, dtype=int), int(iend//99999)+1)
    return ilines[:int(length)].tolist()


def write_line(string, mat, mf, mt, iline):
    return line_pattern.format(string, mat, mf, mt, iline)


def write_eol(lines, mat, mf, mt, istart=1):
    """
    Add end-of-line flags MAT, MF, MT and line number to list of strings.

    Returns
    -------
    `str`
        A string that is the sum of the list of strings, i.e. including
        eol falgs, separated by the newline character `\n`.

    Warns
    -----
    This function does not check if the strings are in ENDF-6 format
    or longer than 66 characters.
    """
    ilines = line_numbers(len(lines))
    return [write_line(string, mat, mf, mt, iline) for string, iline in zip(lines, ilines)]
