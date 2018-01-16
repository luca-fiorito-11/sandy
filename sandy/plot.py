# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:34:12 2017

@author: lfiorito
"""
import matplotlib.pyplot as plt
import logging



def convert_to_stepf(x, y):
    """
    Convert tabulated array into step function for plotting.

    Inputs:
        - x :
            (1d array) x-axis array
        - y :
            (1d array) y-axis array, must have dimensions equal to `len(x)` or 
            `len(x)-1`

    Outputs:
        - X :
            (1d array) step representation of input array x
        - Y :
            (1d array) step representation of input array y
    """
    import numpy as np
    LX = x.size
    LY = LX - 1
    y = y[:LY]
    IX = 2*np.ones(LX, dtype=int)
    IX[0] = IX[-1] = 1
    IY = 2*np.ones(LY, dtype=int)
    X = np.repeat(x, IX, axis=0)
    Y = np.repeat(y, IY, axis=0)
    return X, Y



def save(path, fig, ext='png', close=True, verbose=True):
    """
    Save a figure from pyplot.
    
    Inputs :
        - path :
            (string) the path (and filename, without the extension) to save 
            the figure to
        - ext :
            (string) the file extension.
            This must be supported by the active  matplotlib backend 
            (see matplotlib.backends module).
            Most backends support `png`, `pdf`, `ps`, `eps`, and `svg`
        - close :
            (boolean) whether to close the figure after saving.
            If you want to save the figure multiple times (e.g., to multiple 
            formats), you should NOT close it in between saves or you will 
            have to re-plot it
        - verbose :
            (boolean) whether to print information about when and where the 
            image has been saved
    """
    import os
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "{}.{}".format(os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
    # The final path to save to
    savepath = os.path.join(directory, filename)
    if verbose:
        logging.info("SAVING FIGURE TO '{}'...".format(savepath))
    fig.savefig(savepath, bbox_inches=0)
    if close:
        plt.close(fig)



def plot1d(ax, xx, yy, step=False, **kwargs):
    r"""
    Add plot to `Axes` object in input.
    
    Inputs:
        - ax :
            (`Axes` object)
        - xx :
            (1d-array) tabulated x-values
        - yy :
            (1d-array) tabulated y-values
        - step :
            (boolean) convert tabulated function to step function
    """
    if step:
        XX, YY = convert_to_stepf(xx, yy)
    else:
        XX, YY = xx, yy
    ax.plot(XX, YY, **kwargs)
