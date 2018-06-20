# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 09:12:53 2018

@author: lucaf
"""

import time

def TimeDecorator(foo):
    """
    Output the time a function takes
    to execute.
    """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        out = foo(*args, **kwargs)
        t2 = time.time()
        print("Time to run function {}: {} sec".format(foo, t2-t1))
        return out
    return wrapper
        