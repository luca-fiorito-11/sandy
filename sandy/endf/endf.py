# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 15:01:43 2017

@author: lfiorito
"""

import sys
import logging
import numpy as np
from sandy.files import File


za = 0.
awr = 0.


class CONT(tuple):
    r"""
    The smallest possible record in `ENDF-6` format is a control (`CONT`) 
    record.
    For convenience, a `CONT` record is denoted by::

        [MAT,MF,MT/C1,C2,L1,L2,N1,N2] CONT

    The `CONT` record can be read with the following `FORTRAN` statements::

        READ(LIB,10)C1,C2,L1,L2,N1,N2,MAT,MF,MT,NS
        10 FORMAT(2E11.0,4I11,I4,I2,I3,I5)

    The actual parameters stored in the six fields `C1`, `C2`, `L1`, `L2`, `N1`, 
    and `N2` will depend on the application for the `CONT` record.
    """
    
    def __new__ (cls, string):
        xx = CONT.read_cont(string)
        return tuple.__new__(cls, xx)
    
    @staticmethod
    def read_cont(string):
        """
        Read ``ENDF-6`` ``CONT`` record in formatted fortran.
        
        Outputs:
            - :``out``: :
                (tuple) content of ``CONT`` record
        
        Found error in:
            - n-17-Cl-035.jeff32
            - n-3-Li-007.jeff32
            - n-63-Eu-152.jeff32
            - n-63-Eu-153.jeff32
            - n-64-Gd-155.jeff32
            - n-77-Ir-193.jeff32
            - n-90-Th-229.jeff32
            - n-94-Pu-238.jeff32
            - n-94-Pu-241.jeff32
            - n-94-Pu-242.jeff32
            - n-97-Bk-250.jeff32
            - n-98-Cf-254.jeff32
        """
        import sandy.rw_fortran as rwf
        c1 = np.array(0.0, dtype='d')
        c2 = np.array(0.0, dtype='d')
        l1 = np.array(0, dtype='int')
        l2 = np.array(0, dtype='int')
        n1 = np.array(0, dtype='int')
        n2 = np.array(0, dtype='int')
        mat = np.array(0, dtype='int')
        mf = np.array(0, dtype='int')
        mt = np.array(0, dtype='int')
        ns = np.array(0, dtype='int')
        io_status = np.array(0, dtype='int')
        rwf.rcont(string,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns,io_status)
        if io_status != 0:
            logging.error("ERROR : line '{}' is not in CONT format".format(string))
            sys.exit()
        out = (float(c1), float(c2), int(l1), int(l2), int(n1), int(n2))
        return out
    
    @property
    def c1(self):
        r"""
        First element of string (float).
        """
        return self[0]
    
    @property
    def c2(self):
        r"""
        Second element of string (float).
        """
        return self[1]

    @property
    def l1(self):
        r"""
        Third element of string (integer).
        """
        return self[2]

    @property
    def l2(self):
        r"""
        Fiurth element of string (integer).
        """
        return self[3]

    @property
    def n1(self):
        r"""
        Fifth element of string (integer).
        """
        return self[4]

    @property
    def n2(self):
        r"""
        Sixth element of string (integer).
        """
        return self[5]



class Line(str):
    r"""
    ``ASCII`` textline of the ``ENDF-6`` file.
    """
    
    def __new__(cls, line):
        if len(line) < 80:
            raise NotImplementedError("line '{}' is not 80 charachters long".format(line))
        return str.__new__(cls, line[:80])
    
    @property
    def mat(self):
        """
        Material type (``MAT``).
        """
        return int(self[66:70])
        
    @property
    def mf(self):
        """
        File type (``MF``).
        """
        return int(self[70:72])

    @property
    def mt(self):
        """
        Reaction type (`MT`).
        """
        return int(self[72:75])

    @property
    def ns(self):
        """
        Section number.
        """
        return int(self[75:])
    
    @property
    def CONT(self):
        """
        `ENDF-6` `CONT` record.
        """
        return CONT(self)
    
    @property
    def HEAD(self):
        """
        `ENDF-6` `HEAD` record.
        """
        return self[:66]
        
        


class MF(dict):
    r"""
    Default object to store a `MF` section when the section processing is not
    requested.
    """
    
    def __repr__(self):
        return "<MF{0.mf}>".format(self)
    
    def __init__(self, endf):
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        for mt in sorted(endf.INFO.records[self.mf]):
            if mt != 451:
                self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend
    
    @property
    def prefix(self):
        return "MF{:<2} : ".format(self.mf)
    


class MT:
    r"""
    Default object to store a `MT` section when the section processing is not
    requested.
    """

    def __repr__(self):
        r"""
        `endf.MT` is inherited by the other `MT` sections just to prevent 
        from rewriting magic method `__repr__`.
        """
        return "<MF{0.mf} MT{0.mt}>".format(self)
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.text = self.read_text(endf)
    
    @property
    def prefix(self):
        return "MF{:<2} MT{:<3} : ".format(self.mf, self.mt)
    
    @property
    def za(self):
        r"""
        Floating-point number used to identify materials.
        If `Z` is the charge number and `A` the mass number then `ZA` is 
        computed from
        .. math::
            ZA = (1000.0 \times Z) + A
        """
        return self._za
    
    @za.setter
    def za(self, za):
        from sandy.endf import za as za_reference
        if za != za_reference:
            logging.warn(self.prefix + "found 'za={:.1f}' instead of 'za={:.1f}'".format(za, za_reference))
        self._za = za

    @property
    def awr(self):
        r"""
        Ratio of the mass of the material to that of the neutron.
        """
        return self._awr
    
    @awr.setter
    def awr(self, awr):
        from sandy.endf import awr as awr_reference
        if awr != awr_reference:
            logging.debug(self.prefix + "found 'awr={:.6e}' instead of 'awr={:.6e}'".format(awr_reference, awr))
        self._awr = awr
    
    def read_text(self, endf, getback=False):
        iline = endf.i
        nlines = endf.INFO.records[self.mf][self.mt]
        text = []
        for i in range(nlines):
            line = endf.next().HEAD
            text.append(line)
        line = endf.line
        if line.mat != self.mat or line.mf != self.mf or line.mt != 0:
            logging.error(self.prefix + "Section is not compatible with the number of lines reportes in MF1 MT451: nlines={}".format(nlines))
            sys.exit()
        if getback:
            endf.move_to_line(iline)
        return text
        

    def write(self, **kwargs):
        r"""
        Write the text of the `MT` section as a list of strings.
        Each string corresponds to a `ENDF-6` line, where the control values 
        are updated and the EOL symbol `\n` is added.
        
        .. Important::
            Some objects that inherit ``endf.MF`` redefine the ``write`` 
            method, which could accept extra input arguments.
            ``kwargs`` catches the extra arguments and maintain consistency 
            in the inheritance.
        
        Outputs:
            - text :
                (list of strings) text of the `MT` section
        """
        text = ENDF.add_control(self.text, self.mat, self.mf, self.mt)
        return text



class ENDF(File, dict):
    """
    ``ENDF-6`` control file.
    """
    
    def __repr__(self):
        return "<ENDF-6 file>"

    @classmethod
    def read_pendf(cls, file):
        """
        Given a file in ``ENDF-6`` format, run ``RECONR`` of ``NJOY`` and 
        extract a ``PENDF`` file.
        
        Inputs:
            - :``file``: :
                (``ENDF`` instance) ``ENDF-6`` file
        """
        from sandy.njoy import FileNJOY
        from sandy.sandy_input import outdir, options
        from os.path import join
        from os import unlink
        import shutil
        logging.info("ENDF : Convert ENDF-6 file '{}' to PENDF".format(file.filename))
        fnjoy = FileNJOY('input_njoy')
        # Write NJOY file
        fnjoy.copy_to_tape(file.filename, 20, dst=outdir)
        fnjoy.reconr(20, 21, file.mat, **options)
        fnjoy.stop()
        fnjoy.run(cwd=outdir)
        filependf = join(outdir, file.name + '.pendf')
        shutil.move(join(outdir,'tape21'), filependf)
        # Remove NJOY junk outputs
        unlink(join(outdir,'tape20'))
        unlink(join(outdir,'output'))
        return cls.read_endf(filependf)
    
    @classmethod
    def read_endf(cls, filename):
        """
        Initialize `endf.ENDF` instance starting from a filename.
        
        Inputs:
            - :``filename``: :
                (string) name of the ``ENDF-6`` file
        """
        inst = cls(filename)
        inst.process()
        return inst
 
    @property
    def prefix(self):
        return "ENDF : "
    
    @property
    def title(self):
        r"""
        Title of the ``ENDF-6`` file.

        ..Important::
            Do not create a `endf.Line` object as the length of the title is 
            sometimes :math:`\\leq 80` characters 
        """
        return self.text[0][:66]
    
    @property
    def mat(self):
        r"""
        Material number `MAT`.
        Read it from the second line (after the title) of the `ENDF-6` file.
        """
        return self.get_line(1).mat
    
    @property
    def n(self):
        r"""
        Number of lines of the `ENDF-6` file
        """
        return len(self.text)

    @property
    def line(self):
        """
        Return current line.
        """
        return Line(self.text[self.i])

    def __init__(self, filename):
        logging.debug(self.prefix + "Read file '{}'".format(filename))
        super().__init__(filename)
        self.text = self.read()
        self.i = 0

    def process(self):
        from sandy.sandy_input import options
        from sandy.mf1 import Info
        from sandy.mf1 import MF as MF1
        from sandy.mf2 import MF as MF2
        from sandy.mf3 import MF as MF3
        from sandy.mf4 import MF as MF4
        from sandy.mf5 import MF as MF5
        from sandy.mf32 import MF as MF32
        from sandy.mf33 import MF as MF33
        from sandy.mf34 import MF as MF34
        from sandy.mf35 import MF as MF35
        self.move_to_line(1) # skip title
        self.INFO = Info(self)
        for mf in sorted(self.INFO.records):
            if mf == 1 and mf in options['mf']:
                self[mf] = MF1(self)
            elif mf == 2 and mf in options['mf'] and self.INFO.lrp != 2:
                self[mf] = MF2(self)
            elif mf == 3 and mf in options['mf']:
                self[mf] = MF3(self)
            elif mf == 4 and mf in options['mf']:
                self[mf] = MF4(self)
            elif mf == 5 and mf in options['mf']:
                self[mf] = MF5(self)
            elif mf == 31 and mf in options['mf']:
                self[mf] = MF33(self)
            elif mf == 32 and mf in options['mf']:
                self[mf] = MF32(self)
            elif mf == 33 and mf in options['mf']:
                self[mf] = MF33(self)
            elif mf == 34 and mf in options['mf']:
                self[mf] = MF34(self)
            elif mf == 35 and mf in options['mf']:
                self[mf] = MF35(self)
            else:
                self[mf] = MF(self)

    def __iter__(self):
        """
        Make the file an iterator.
        """
        return self
    
    def next(self):
        """
        Yield a new line of the `ENDF-6` file.
        """
        if self.i < self.n:
            line = self.line
            self.i += 1
            return line
        else:
            raise NotImplementedError("Reached EOF of the `ENDF-6` file")
    
    def back(self, steps=1):
        """
        Go back a number of lines in the `ENDF-6` file.
        
        Inputs:
            - steps :
                (scalar integer) number of lines to jump backword, default is 
                one.
        """
        self.i -= steps
        if self.i <= 0:
            raise NotImplementedError("Reached BOF of the `ENDF-6` file")
        
    def move_to_line(self, idx):
        """
        Given an index in input, move pointer to the corresponding line.
        
        Inputs:
            - idx :
                (scalar integer) line position
        """
        i = int(idx)
        if i < 0:
            raise NotImplementedError("Cannot move to negative lines in `ENDF-6` file")
        if i >= self.n:
            raise NotImplementedError("Trying to move beyond EOF of the `ENDF-6` file")
        self.i = i
    
    
    def get_line(self, idx=None):
        """
        Return line corresponing to index position.
        
        Inputs:
            - idx :
                (scalar integer) line position, default is the pointer
        """
        if idx is None:
            return self.line
        i = int(idx)
        if i < 0:
            raise NotImplementedError("Cannot get negative lines in `ENDF-6` file")
        if i >= self.n:
            raise NotImplementedError("Trying to get line beyond EOF of the `ENDF-6` file")
        return Line(self.text[i])
            
    
    def read_control(self):
        """
        Read `ENDF-6` control record, i.e. the numbers in each ENDF-6 file line 
        after column 66.
        
        Outputs:
            - mat :
                (scalar integer) material identifier
            - mf :
                (scalar integer) file identifier
            - mt :
                (scalar integer) section identifier
            - ns :
                (scalar integer) line number
        """
        line = self.line
        mat = line.mat
        mf = line.mf
        mt = line.mt
        ns = line.ns
        return mat, mf, mt, ns

    
    def read_text(self):
        """
        Read `ENDF-6` `HEAD` record (used to give a comment or title).
        
        Outputs:
            - head : 
                (string) first 66 characters of the line
        """
        head = self.next().HEAD
        return head

    
    def read_cont(self, delete=[]):
        """
        Read `ENDF-6` `CONT` record in formatted fortran.
        
        Inputs:
            - delete :
                (list of integers) indexes of the elements that must be 
                deleted in `CONT` record (the default is to keep all)
        Outputs:
            - cont :
                (tuple) `endf.CONT` object
        """
        cont = list(self.next().CONT)
        for i in sorted(delete, reverse=True):
            del cont[i]
        cont = tuple(cont)
        if len(cont) == 1:
            cont = cont[0]
        return cont


    def read_list(self):
        """
        Read ENDF-6 **list** record.
        
        Outputs:
            - c1 : 
                (float) first element of string
            - c2 : 
                (float) second element of string
            - l1 : 
                (int) third element of string
            - l2 : 
                (int) fourth element of string
            - npl : 
                (int) number of elements in the list
            - n2 : 
                (int) sixth element of string
            - b :
                (list) series of numbers
        """
        import sandy.rw_fortran as rwf
        c1, c2, l1, l2, npl, n2 = self.read_cont()
        count = 0
        b = []
        while count < npl:
            iline = self.i
            string = self.next()
            array = np.zeros(6, dtype=float)
            io_status = np.array(0, dtype='int')
            rwf.rlist(string, io_status, array, 6)
            if io_status != 0:
                logging.error("ERROR! Cannot read line {} of file '{}'".format(iline, self.filename))
                sys.exit()
            b.extend(array)
            count += 6            # In fortran fields are read 6 by 6
        return c1, c2, l1, l2, npl, n2, b[:npl]


    def read_tab2(self):
        """
        Read ENDF-6 **tab2** record.
        This record is used to control the tabulation of a 2-dimensional 
        function :math:`y(x,z)`.
        It specifies how many values of z are to be given and how to
        interpolate between the successive values of :math:`z`.
        Tabulated values of :math:`y_l(x)` at each value of :math:`z_l` 
        are given in ENDF-6 **tab1** or **list** records following the 
        ENDF-6 **tab2** record, with the appropriate value of :math:`z` 
        in the field designated as *c2*.
        
        Outputs:
            - c1 : 
                (float) first element of string
            - c2 : 
                (float) second element of string
            - l1 : 
                (int) third element of string
            - l2 : 
                (int) fourth element of string
            - nr :
                (int) number of interpolation ranges
            - nz :
                (int) number of :math:`z` to be given
            - nbt :
                (list of int) upper limits of the interpolation ranges
            - interp :
                (list of int) interpolation functions
        """
        import sandy.rw_fortran as rwf
        c1, c2, l1, l2, nr, nz = self.read_cont()
        count = 0
        tab = []
        while count < nr*2:
            iline = self.i
            string = self.next()
            array = np.zeros(6, dtype=float)
            io_status = np.array(0, dtype='int')
            rwf.rlist(string, io_status, array, 6)
            if io_status != 0:
                logging.error("ERROR! Cannot read line {} of file '{}'".format(iline, self.filename))
                sys.exit()
            tab.extend(list(map(int, array)))
            count += 6
        nbt = tab[::2][:nr]
        interp = tab[1::2][:nr]
        return c1, c2, l1, l2, nr, nz, nbt, interp
    

    def read_tab1(self):
        """
        Read ENDF-6 **tab1** record.
        These records are used for one-dimensional tabulated functions such as 
        :math:`y(x)`.
        The data needed to specify a one-dimensional tabulated function are 
        the interpolation tables *nbt* and *int* for each of the *nr* ranges, 
        and the *np* tabulated pairs of :math:`x(n)` and :math:`y(n)`.
        
        Outputs:
            - c1 : 
                (float) first element of string
            - c2 : 
                (float) second element of string
            - l1 : 
                (int) third element of string
            - l2 : 
                (int) fourth element of string
            - nr :
                (int) number of interpolation ranges
            - nz :
                (int) number of :math:`z` to be given
            - nbt :
                (list of int) upper limits of the interpolation ranges
            - interp :
                (list of int) interpolation functions
            - x, y : 
                tabulated pairs of :math:`x(n)` and :math:`y(n)`
        """
        import sandy.rw_fortran as rwf
        c1, c2, l1, l2, nr, nz, nbt, interp = self.read_tab2()
        count = 0
        tab = []
        while count < nz*2:
            iline = self.i
            string = self.next()
            array = np.zeros(6, dtype=float)
            io_status = np.array(0, dtype='int')
            rwf.rlist(string, io_status, array, 6)
            if io_status != 0:
                logging.error("ERROR! Cannot read line {} of file '{}'".format(iline, self.filename))
                sys.exit()
            tab.extend(array)
            count += 6            # In fortran fields are read 6 by 6
        x = tab[::2][:nz]
        y = tab[1::2][:nz]
        return c1, c2, l1, l2, nr, nz, nbt, interp, x, y

    def read_compact_cov(self):
        import sandy.rw_fortran as rwf
        # Read the CONT record:
        # NNN is the dimension of CORR(NNN,NNN),
        # NM is the number of lines to follow in the file
        # NDIGIT is the number of digits for the covariance matrix
        c1, c2, ndigit, nnn, nm, nx = self.read_cont()
        # Preset the correlation matrix to zero
        # Preset the diagonal to one       
        matrix = np.eye(nnn, dtype='float')
        for m in range(nm):   #        Read the INTG record
            iline = self.i
            string = self.next()
            kij = np.zeros(18, dtype='float')
            ii = np.array(0, dtype='int')
            jj = np.array(0, dtype='int')
            nrow = np.array(0,dtype='int')
            io_status = np.array(0, dtype='int')
            rwf.lcomp(ndigit, string, io_status, ii, jj, nrow, kij, 18)
            if io_status != 0:
                logging.error("ERROR! Cannot read line {} of file '{}'".format(iline, self.filename))
                sys.exit()
            ii = int(ii) - 1
            jj = int(jj) - 1
            nrow = int(nrow)
            kij.astype(int)
            jp = jj - 1
            Factor = 10**(ndigit)
            for n in range(nrow):
                jp += 1
                if jp >= ii:
                    break
                if kij[n] != 0:
                    if kij[n] > 0:
                        cij = (kij[n] + 0.5)/Factor
                    else:
                        cij = -(-kij[n] + 0.5)/Factor
                    if ii >= nnn:
                        logging.error('COMP2 : Index {} exceeds matrix dimension {}'.format(ii,nnn))
                        sys.exit()
                    if jp >= nnn:
                        logging.error('COMP2 : Index {} exceeds matrix dimension {}'.format(jp,nnn))
                        sys.exit()
                    matrix[ii,jp] = matrix[jp,ii] = cij
        return matrix

    
    @staticmethod
    def write_cont(c1, c2, l1, l2, n1, n2):
        """
        Write ENDF-6 **cont** record.
        
        Outputs:
            - list of string
        """
        import sandy.rw_fortran as rwf
        c1 = np.array(c1, dtype=float)
        c2 = np.array(c2, dtype=float)
        l1 = np.array(l1, dtype=int)
        l2 = np.array(l2, dtype=int)
        n1 = np.array(n1, dtype=int)
        n2 = np.array(n2, dtype=int)
        string = np.array("*"*67)
        rwf.wcont(string, c1, c2, l1, l2, n1, n2)
        string = str(string, 'utf-8')[:66] # byte string coming from f2py must be converted
        return [string]
    
    @staticmethod
    def spaces_to_eol(string):
        if len(string) == 66:
            return string
        else:
            return string + ' '*(66-len(string))
            
    @staticmethod
    def wlist(b, dtype=float):
        import sandy.rw_fortran as rwf
        from sandy.functions import split_by_n
        length = len(b)*11
        string = np.array("*"*(length + 1))
        if dtype is int:
            rwf.wlist_int(string, b, length)
        else:
            rwf.wlist_float(string, b, length)
        string = str(string, 'utf-8')[:length]
        list_of_strings = list(split_by_n(string, 66))
        list_of_strings = [ ENDF.spaces_to_eol(line) 
                            for line in list_of_strings]
        return list_of_strings

    @staticmethod
    def write_list(c1, c2, l1, l2, n2, b):
        """
        Write ENDF-6 **list** record.
        
        Outputs :
            list of strings.
        """
        strings = ENDF.write_cont(c1, c2, l1, l2, len(b), n2)
        list_of_strings = ENDF.wlist(b)
        strings.extend(list_of_strings)
        return strings
    
    @staticmethod
    def write_tab2(c1, c2, l1, l2, n1, n2, nbt, intr):
        """
        Write ENDF-6 **tab2** record.

        Outputs:
            list of strings
        """
        strings = ENDF.write_cont(c1, c2, l1, l2, n1, n2)
        b = [ x for item in zip(nbt, intr) for x in item ]
        list_of_strings = ENDF.wlist(b, int)
        strings.extend(list_of_strings)
        return strings

    @staticmethod
    def write_tab1(c1, c2, l1, l2, nbt, intr, x, y):
        """
        Write ENDF-6 **tab1** record.
        
        Outputs:
            list of strings
        """
        n1 = len(nbt)
        n2 = len(x)
        strings = ENDF.write_cont(c1, c2, l1, l2, n1, n2)
        b = [ z for item in zip(nbt, intr) for z in item ]
        list_of_strings = ENDF.wlist(b, int)
        strings.extend(list_of_strings)
        b = [ z for item in zip(x, y) for z in item ]
        list_of_strings = ENDF.wlist(b, float)
        strings.extend(list_of_strings)
        return strings

    @staticmethod
    def add_control(text, mat, mf, mt, ns=1):
        """
        Add control values `MAT`, `MF`, `MT` and line number to each line 
        of the input text.
        
        Inputs:
            - mat:
                (scalar integer) material number
            - mf :
                (scalar integer) file number
            - mt :
                (scalar integer) reaction number
            - ns :
                (scalar integer) line number, by default start counting from 1
        
        Outputs:
            - new_text :
                (list of strings) text with control values and EOL symbol `\n`
        """
        new_text = []
        for line in text:
            if ns > 99999:
                ns = 1
            new_line = line[:66] + "{:>4}{:>2}{:>3}{:>5}\n".format(mat, mf, mt, ns)
            new_text.append(new_line)
            ns += 1
        return new_text

    @staticmethod
    def write_send(mat, mf):
        line = ENDF.write_cont(0, 0, 0, 0, 0, 0)
        line = ENDF.add_control(line, mat, mf, 0, ns=99999)
        return line

    @staticmethod
    def write_fend(mat):
        line = ENDF.write_cont(0, 0, 0, 0, 0, 0)
        line = ENDF.add_control(line, mat, 0, 0, ns=0)
        return line

    @staticmethod
    def write_mend():
        line = ENDF.write_cont(0, 0, 0, 0, 0, 0)
        line = ENDF.add_control(line, 0, 0, 0, ns=0)
        return line

    @staticmethod
    def write_tend():
        line = ENDF.write_cont(0, 0, 0, 0, 0, 0)
        line = ENDF.add_control(line, -1, 0, 0, ns=0)
        return line

    def write(self, ismp=None):
        """
        Print file text in a list of stings using the `ENDF-6` format rules.
        
        Inputs:
            - ismp :
                (scalar integer) sample index, default is no samples, write 
                best-estimates

        Outputs:
            - text :
                (list of string) `ENDF-6` file text
        """
        # NOW FILL IN THE TEXT LIST
        text = []
        for mf,mfsec in sorted(self.items()):
#            if mf >= 30:
#                del self[mf]
#                continue
            if len(mfsec) == 0:
                continue
            for mt,mtsec in sorted(mfsec.items()):
                mtsec.text = mtsec.write(ismp=ismp)
                text += mtsec.text
                text += self.write_send(self.mat, mf)
            text += self.write_fend(self.mat)
        text += self.write_mend()
        text += self.write_tend()
        text1451 = self.INFO.write(self)
        text1451 += self.write_send(self.mat, 1)
        if len(self.INFO.records[1]) == 1:
            text1451 += self.write_fend(self.mat)
        text = text1451 + text
        text = self.add_control([self.title], 0, 0, 0, 0) + text
        return text
    
    def write_to_file(self, file, ismp=None):
        """
        Write data to file.
        
        Inputs:
            - :``file``: :
                output file
            - :``ismp``: :
                index of the sample that must be written (default do not 
                write samples, write best estimates)
        """
        text = self.write(ismp=ismp)
        with open(file, 'w', encoding="ascii") as f:
            string = "".join([ line.encode('ascii', 'replace').decode('ascii') for line in text ])
            f.write(string)
