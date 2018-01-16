# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:35:07 2017

@author: lfiorito
"""
import numpy as np
import logging

def triu_matrix(arr, size):
    """
    Given the upper triangular values of a **square symmetric** matrix in 
    an array, return the full matrix.
    
    Inputs:
        - arr :
            (1d array) array with the upper triangular values of the matrix
        - size :
            (int) dimension of the matrix
    
    Outputs:
        - matrix :
            (2d array) reconstructed 2d-array with symmetric matrix
    """
    matrix = np.zeros([size, size])
    indices = np.triu_indices(size)
    matrix[indices] = arr
    matrix += np.triu(matrix, 1).T
    return matrix



def up2down(C):
    """
    Given a covariance matrix in input, copy the upper triangular part to the 
    lower triangular part.
    
    Inputs:
        - C : 
            (2d-array) input covariance matrix
    
    Outputs:
        - C1 : 
            (2d-array) output covariance matrix
    """
    mytype = type(C)
    U = np.triu(C)
    L = np.triu(C, 1).T
    C1 = mytype(U + L)
    return C1



def corr2cov(corr, s):
    mytype = type(corr)
    dim = corr.shape[0]
    S = np.repeat(s, dim).reshape(dim, dim)
    cov = mytype(S.T * (corr * S))
    cov = up2down(cov)
    return cov



class Cov(np.ndarray):
    
    def __new__(cls, arr):
        obj = np.ndarray.__new__(cls, arr.shape, arr.dtype)
        obj[:] = arr[:]
        return obj
    
    @property
    def prefix(self):
        return "COV : "

    @property
    def dim(self):
        """
        Dimension of the covariance.
        """
        length = self.shape[0]
        return length

    @property
    def var(self):
        r"""
        Variance array.
        """
        var = np.diag(np.array(self))
        return var

    @property
    def std(self):
        r"""
        Standard deviation array.
        """
        var = self.var
        if (var < 0).any():
            raise ValueError("Variances must be non-negative")
        std = np.sqrt(var)
        return std
    
    @property
    def nnegvar(self):
        r"""
        Number of negative variances.
        """
        return np.flatnonzero(self.var < 0).size

    @property
    def nzerovar(self):
        r"""
        Number of zero variances.
        """
        return np.flatnonzero(self.var == 0).size
    
    def empty_off_diagonals(self):
        r"""
        Remove off-diagonal elements.
        
        Outputs:
            - :``C``: :
                (``cov.Cov instance``) covariance with empty off-diagonals
        """
        logging.info(self.prefix + "'no_correlations' option is requested, delete off-diagonal terms")
        C = Cov(np.diag(np.diag(self)))
        return C

    def is_symmetric(self):
        r"""
        Check if covariance is symmetric.
        
        If it is nearly symmetric (rtol=1e-5), then we copy the upper 
        triangular part to the lower triangular part and we make it 
        symmetric.
        
        Outputs:
            - :``check``: :
                (boolean) ``True`` if matrix is symmetric, else ``False``
        """
        check = True
        if not (self.T == self).all():
            check = False
            if np.isclose(self.T, self).all():
                check = True
                self[:] = up2down(self)
        return check
    
    def reduce_size(self):
        """
        Reduce matrix dimensions when zeros are found on the diagonal.
        
        Outputs:
            * :``nonzero_idxs``: :
                (1d array) positions of the original diagonal matrix where the 
                coefficients were not zero
            * :``cov_reduced``: :
                (``cov.Cov`` instance) reduced covariance matrix
                
        """
        nonzero_idxs =  np.flatnonzero(np.diag(self))
        cov_reduced = self[nonzero_idxs][:,nonzero_idxs]
        return nonzero_idxs, cov_reduced

    def restore_size(self, nonzero_idxs, cov_reduced):
        """
        Restore original matrix dimensions from a reduced matrix and an array 
        of positions to convert from reduced to original size.
        
        Inputs:
            * :``nonzero_idxs``: :
                (1d array) positions of the original diagonal matrix where the 
                coefficients were not zero
            * :``cov_reduced``: :
                (``cov.Cov`` instance) reduced covariance matrix
        
        Outputs:
            * :``cov``: :
                (``cov.Cov`` instance) reduced covariance matrix increased 
                to given size according to the indexes given in input
        """
        cov = Cov(np.zeros_like(self))
        for i,ni in enumerate(nonzero_idxs):
            cov[ni,nonzero_idxs] = cov_reduced[i]
        return cov

    def sampling(self, nsmp, pdf='normal'):
        r"""
        Extract random samples from the covariance matrix, either using 
        the cholesky or the eigenvalue decomposition.

        Inputs:
            - :``nsmp``: :
                (integer) number of samples
         
        Outputs:
            - :``samples``: :
                (array) random samples
        """
        from sandy.sandy_input import options
        logging.debug(self.prefix + "Covariance matrix dimension is {} X {}".format(*self.shape))
        y = np.random.randn(self.dim, int(nsmp))
        nonzero_idxs, cov_reduced = self.reduce_size()
        nzeros = self.shape[0] - len(nonzero_idxs)
        if nzeros > 0:
            logging.debug(self.prefix + "Found {} zeros on the diagonal, reduce matrix dimension to {} X {}".format(nzeros, *cov_reduced.shape))
        if "no_correlations" in options:
            if options['no_correlations']:
                cov_reduced = cov_reduced.empty_off_diagonals()
        try:
            L_reduced = cov_reduced.cholesky()
            logging.debug(self.prefix + "Cholesky decomposition was successful")
        except np.linalg.linalg.LinAlgError as exc:
            logging.debug(self.prefix + "Cholesky decomposition was not successful, proceed with eigenvalue decomposition")
            L_reduced = cov_reduced.eigendecomp()
        L = self.restore_size(nonzero_idxs, L_reduced)
        samples = np.array(L.dot(y), dtype=float)
        return samples
    
    @property
    def corr(self):
        r"""
        Correlation matrix.
        """
        from sandy.functions import div0
        if not self.is_symmetric():
            raise ValueError("Covariance matrix must be square and symmetric")
        coeff = div0(1, self.std)
        corr = np.multiply(np.multiply(self.T, coeff).T, coeff)
        return corr

    def cholesky(self):
        r"""
        Perform a Cholesky decomposition of the covariance matrix.

        Outputs:
            - :``L``: : 
                (2d array) lower triangular matrix
        """
        from scipy.linalg import cholesky
        L = cholesky(self, lower=True, overwrite_a=False, check_finite=False)
        return L

    def eig(self):
        r"""
        Extract eigenvalues and eigenvectors of the covariance matrix.

        Outputs:
            - :``E``: :
                (1d-array) eigenvalues
            - :``V``: :
                (2d-array) eigenvectors
        """
        from scipy.linalg import eig
        E, V = eig(self)
        E, V = E.real, V.real
        return E, V
    
    def eigendecomp(self):
        r"""
        Perform an eigenvalue decomposition of the covariance matrix.

        Outputs:
            - :``L``: :
                (2d-array) lower triangular matrix
        """
        from scipy.linalg import qr
        E, V = self.eig()
        negative_eig = np.extract(E < 0, E)    # extract negative eigenvalues
        if len(negative_eig) != 0:
            largest_negative = max(abs(negative_eig))
            logging.debug(self.prefix + '{} negative eigenvalues were found and replaced with zero'.format(negative_eig.size))
            pos = sorted(abs(E),reverse=True).index(largest_negative) + 1
            logging.debug(self.prefix + 'Largest negative eigenvalue ranks {}/{}'.format(pos, E.size))
            logging.debug(self.prefix + 'eig(-)/eig_max = {}%'.format(largest_negative/max(abs(E))*100.))
        E[E<=0] = 0
        Esqrt = np.diag(np.sqrt(E))
        M = V.dot(Esqrt)
        Q,R = qr(M.T)
        L = R.T
        logging.debug(self.prefix + "Eigenvalue decomposition was successful")
        return L
    
    def plot(self):
        r"""
        Plot covariance matrix as a pseudocolor plot of a 2-D array.
        The colorbar is also added to the figure.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        pcm = ax.matshow(self.corr, vmin=-1, vmax=1, cmap='bwr', aspect='auto')
        # Resize the plot to make space for the colorbar
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, 0.7, box.height])
        # set labels
        ax.set_title('evaluated correlation matrix')
        ax.set_xlabel('energy (eV)')
        ax.set_ylabel('energy (eV)')
        # Plot the colorbar in desired position 
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8]) 
        plt.colorbar(pcm, cax=cbaxes)
        fig.show()
    
    def dump(self, fname):
        np.savetxt(fname, self, fmt='%.5e')



class Ucov(dict):
    r"""
    
    *data* must be structured as a dictionary::
        
        data = { (key1,key2) : Cov33, ... }
    
    Hence, each dictionary keyword must be a tuple with two elements.
    Furthermore, each key can also be a tuple or else (not ndarray), 
    as long as *key1* and *key2* remain symmetric.
        
    Extract pointers for the union covariance matrix.
    
    Example,
    
    |              l=0   l=0   l=0   l=1   l=1
    |              e=1   e=2   e=3   e=2   e=3
    | l=0 e=1       x     x     x     x     x
    | l=0 e=2       x     x     x     x     x
    | l=0 e=3       x     x     x     x     x 
    | l=1 e=2       x     x     x     x     x
    | l=1 e=3       x     x     x     x     x
    
    Then the pointers are `[0,0,0,1,1]`.
    
    Inputs:
        - mt1 : 
            (int) "other" reaction type
    
    Ouputs:
        - pointers :
            1d array of pointers
    """
    def __setitem__(self, key, item):
        """
        Filter the items to set in a ``cov.Ucov`` dictionary according to 
        the following rules:
            * raise error if key is not a tuple of two elements
            * raise error if item is not of type ``mf33.Cov33``
        """
        from sandy.mf33 import Cov33
        if len(key) != 2:
            raise NotImplementedError("'Ucov' keys must be a tuple with size=2")
        if isinstance(item, Cov33):
            raise NotImplementedError("Item for 'key={}' must be of type 'Cov33'".format(key))
        super().__setitem__(key, item)

    def get_ucov(self):
        from collections import namedtuple
        Pointer = namedtuple('Pointer', ('beg', 'end'))
        indices = []
        for key in sorted(self):
            i,j = key
            if i != j:
                continue
            indices += [i]*len(self[i,j].xx)
        indices = np.array(indices)
        cov = Cov(np.zeros((len(indices),)*2))
        pntrs = {}
        for i,j in sorted(self):
            # this is for tuple indices
            mask = np.array([(ind == i).all() for ind in indices])
            beg = (mask).nonzero()[0][0]
            end = (mask).nonzero()[0][-1]
            pntrs[i] = Pointer(beg, end)
        for i,j in sorted(self):
            xi = self[i,i].xx
            xj = self[j,j].xx
            C = self[i,j].reinit(xi, xj)
            cov[pntrs[i].beg:pntrs[i].end+1, pntrs[j].beg:pntrs[j].end+1] = C
            cov[pntrs[j].beg:pntrs[j].end+1, pntrs[i].beg:pntrs[i].end+1] = C.T
        return pntrs, cov

    def sampling(self, nsmp):
        from sandy.sandy_input import options
        from sandy.records import Samples
        pntrs, cov = self.get_ucov()
        if "no_correlations" in options:
            if options['no_correlations']:
                cov = cov.empty_off_diagonals()
        samples = cov.sampling(nsmp)
        DICT = {}
        for k,p in pntrs.items():
            xx = self[k,k].xx
            DICT[k] = Samples(xx, samples[p.beg:p.end+1])
        return DICT

    def plot(self):
        """
        Plot the correlation matrix of the union covariance as a pseudocolor 
        plot of a 2-D array.
        The colorbar is also added in the figure.
        """
        import matplotlib.pyplot as plt
        pntrs, cov = self.get_ucov()
        fig, ax = plt.subplots()
        corr = cov.corr()
        pcm = ax.matshow(corr, vmin=-1, vmax=1, cmap='bwr', aspect='auto')
        
        # now, add spacers to separate the different keys
        xtick_name = []
        ytick_name = []
        tick_pos = []
        
        for mt,p in sorted(pntrs.items()):
            ax.axhline(y=p.beg-0.5, linewidth=1, color='k')
            ax.axvline(x=p.beg-0.5, linewidth=1, color='k')
            delta = p.end - p.beg
            tick_pos.append(p.beg + delta/2)
            ytick_name.append("{}".format(mt))
            xtick_name.append("{}".format(mt))

        plt.xticks(tick_pos, xtick_name)
        plt.yticks(tick_pos, xtick_name)
        cbar = plt.colorbar(pcm)
        cbar.ax.tick_params()
        fig.show()



