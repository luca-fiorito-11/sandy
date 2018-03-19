# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 23:41:26 2018

@author: lucaf
"""
import multiprocessing as mp

def cube(x):
    return x**3

if __name__ == '__main__':  
    pool = mp.Pool(processes=4)
    results = [pool.apply(cube, args=(x,)) for x in range(1,7)]
    print(results)