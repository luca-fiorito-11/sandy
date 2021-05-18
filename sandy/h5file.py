import h5py
import os

import sandy


__author__ = "Luca Fiorito"

__all__ = [
        "H5File",
        ]


class H5File():

    def __init__(self, file):
        self.file = file

    @property
    def file(self):
        """
        HDF5 filename with absolute or relative path.

        Returns
        -------
        `str`
            filename
        """
        return self._file

    @file.setter
    def file(self, file):
        self._file = file

    @property
    def data(self):
        """
        HDF5 file content accessible via `h5py.File` interface.

        Returns
        -------
        `h5py.File`
            content of HDF5 file

        Raises
        ------
        `aleph.Error`
            if file is not open

        Notes
        -----
        .. important:: Attribute `data` is allocated any time method `open`
                       is called, and deallocated when method `close` is
                       called.
        """
        if not hasattr(self, "_data"):
            raise sandy.Error("HDF5 file must be opened first")
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

    def open(self, mode="r", **kwargs):
        """
        Open HDF5 file using `h5py.File`.

        Parameters
        ----------
        mode : `str`, optional, default is `'r'`
            opening mode, default is *read only*
        **kwargs : `dict`, optional
            keywrod arguments to pass to `h5py.File`

        Notes
        -----
        .. note:: This method works **inplace**: the `h5py.File` instance
                  is stored as attribute `data`.
        """
        self.data = h5py.File(self.file, mode=mode, **kwargs)

    def close(self, **kwargs):
        """
        Close HDF5 file.

        Parameters
        ----------
        **kwargs : `dict`, optional
            keywrod arguments to pass to `h5py.File.close`

        Notes
        -----
        .. note:: Any time this method is called, attribute `data` is
                  deallocated.
        """
        if hasattr(self, "_data"):
            self.data.close(**kwargs)
            del self._data

    def exists(self):
        """
        Check if file already exists.

        Returns
        -------
        `bool`
            `True` if file exists, else `False`
        """
        return os.path.exists(self.file)
