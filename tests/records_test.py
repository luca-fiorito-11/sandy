import pytest

import sandy


@pytest.mark.parametrize(
        "lines, expected",
        [
                ([' 3.50000000-4.00000000          6          7          8     -64649'], (3.5, -4, 6, 7, 8, -64649)),
                ]
        )
def test_read_cont(lines, expected):
    ipos = 0
    C, ipos = sandy.read_cont(lines, ipos)
    assert ipos == 1
    assert C.C1 == 3.5
    assert C.C2 == -4
    assert C.L1 == 6
    assert C.L2 == 7
    assert C.N1 == 8
    assert C.N2 == -64649



@pytest.mark.parametrize(
        "value, expected",
        [
                (1.123456789e11,' 1.12346+11'),
                (1.123456789e10,' 1.12346+10'),
                (1.123456789e9, ' 1123456789'),
                (1.123456789e8, '  112345679'),
                (1.123456789e7, ' 11234567.9'),
                (1.123456789e6, ' 1123456.79'),
                (1.123456789e5, ' 112345.679'),
                (1.123456789e4, ' 11234.5679'),
                (1.123456789e3, ' 1123.45679'),
                (1.123456789e2, ' 112.345679'),
                (1.123456789e1, ' 11.2345679'),
                (1.123456789e0, ' 1.12345679'),
                (0  , ' 0.00000000'),
                (1.123456789e-1,' 1.123457-1'),
                (1.123456789e-2,' 1.123457-2'),
                (1.123456789e-3,' 1.123457-3'),
                (1.123456789e-4,' 1.123457-4'),
                (1.123456789e-5,' 1.123457-5'),
                (1.123456789e-6,' 1.123457-6'),
                (1.123456789e-7,' 1.123457-7'),
                (1.123456789e-8,' 1.123457-8'),
                (1.123456789e-9,' 1.123457-9'),
                (1.123456789e-10, ' 1.12346-10'),
                (1.123456789e-11, ' 1.12346-11'),
                (-1e11,'-1.00000+11'),
                (-1e10,'-1.00000+10'),
                (-1e9, '-1000000000'),
                (-1e8, ' -100000000'),
                (-1e7, '-10000000.0'),
                (-1e6, '-1000000.00'),
                (-1e5, '-100000.000'),
                (-1e4, '-10000.0000'),
                (-1e3, '-1000.00000'),
                (-1e2, '-100.000000'),
                (-1e1, '-10.0000000'),
                (-1e0, '-1.00000000'),
                (-1e-1, '-1.000000-1'),
                (-1e-2, '-1.000000-2'),
                (-1e-3, '-1.000000-3'),
                (-1e-4, '-1.000000-4'),
                (-1e-5, '-1.000000-5'),
                (-1e-6, '-1.000000-6'),
                (-1e-7, '-1.000000-7'),
                (-1e-8, '-1.000000-8'),
                (-1e-9, '-1.000000-9'),
                (-1e-10,'-1.00000-10'),
                (-1e-11,'-1.00000-11'),
                ]
        )
def test_write_float(value, expected):
    out = sandy.write_float(value)
    assert out == expected



@pytest.mark.parametrize(
        "values, expected",
        [
                ((3.5, -4, 6, 7, 8, -64649), ' 3.50000000-4.00000000          6          7          8     -64649'),
                ((3.5, -4, 6, 7.5, 8, -64649), Exception),
                ((3.5, -4, 6, 7, 8, "aa"), Exception),
                ]
        )
def test_write_cont(values, expected):
    if expected is Exception:
        with pytest.raises(Exception):
            sandy.write_cont(*values)
    else:
        lines = sandy.write_cont(*values)
        assert lines[0] == expected



@pytest.mark.parametrize(
        "lst, expected",
        [
                (list(range(10)), ['          0          1          2          3          4          5',
                                   '          6          7          8          9                      ']),
                ([10], ['         10                                                       ']),
                (10, Exception),
                ([3.5], Exception),
                (["aaa"], Exception),
                ]
        )
def test_write_integer_list(lst, expected):
    if expected is Exception:
        with pytest.raises(Exception):
            sandy.write_integer_list(lst)
    else:
        lines = sandy.write_integer_list(lst)
        assert lines == expected



@pytest.mark.parametrize(
        "lst, expected",
        [
                (list(range(10)), [' 0.00000000 1.00000000 2.00000000 3.00000000 4.00000000 5.00000000',
                                   ' 6.00000000 7.00000000 8.00000000 9.00000000                      ']),
                ([10.5], [' 10.5000000                                                       ']),
                (10.5, Exception),
                (["aaa"], Exception),
                ]
        )
def test_write_float_list(lst, expected):
    if expected is Exception:
        with pytest.raises(Exception):
            sandy.write_float_list(lst)
    else:
        lines = sandy.write_float_list(lst)
        assert lines == expected



@pytest.mark.parametrize(
        "length, istart, expected",
        [
                (10, 2, [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
                (10, 99997, [99997, 99998, 99999, 1, 2, 3, 4, 5, 6, 7]),
                (3, -2, [-2, -1, 0]),
                ]
        )
def test_line_numbers(length, istart, expected):
    lines = sandy.line_numbers(length, istart=istart)
    assert lines == expected



@pytest.mark.parametrize(
        "lines, mat, mf, mt, expected",
        [
                (["aaa","bbb"], 9237, 3, 2, ['aaa                                                               9237 3  2    1',
                                             'bbb                                                               9237 3  2    2']),
                (["*"*100], 9237, 3, 2, ["*"*100+"9237 3  2    1"]),
                (["aaa","bbb"], 9237, 3.5, 2, Exception),
                (["aaa","bbb"], 9237.1, 3, 2, Exception),
                (["aaa","bbb"], 9237, 3, 2.0, Exception),
                ]
        )
def test_write_eol(lines, mat, mf, mt, expected):
    if expected is Exception:
        with pytest.raises(Exception):
            sandy.write_eol(lines, mat, mf, mt)
    else:
        string = sandy.write_eol(lines, mat, mf, mt)
        assert string == expected