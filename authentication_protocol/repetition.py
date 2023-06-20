# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 21:34:43 2020

@author: chrys
"""


#!/usr/bin/env python3

# Chrys 2020-12-15

#standard libary
from math import *
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
#mylibrary
import myfunc as myf

if __name__ == "__main__":
    #it give the probability of generating an identical key after repetition coding.
    #OLD See version1 for improved graph
    
    #variables
    l = 100 #length of correlated sequences
    M = [1,2,3,4,5] #blocksize of repetition coding
    A = np.linspace(0,0.5,num=50) #probability of mismatch
    
    
    ExpErrorsAfterRepet = [int(l/2)* a**m for m in M  for a in A]
    Ma = np.resize(ExpErrorsAfterRepet, (len(M), len(A))) # different row, different m


    #graph heatmap
    heat = plt.pcolormesh(A, M, Ma)
    plt.colorbar(heat)
    #plt.plot(A,ExpErrorsAfterRepet)
    #plt.legend(["Expected errors after repetition"])
    #plt.title("Pr(bit mismatch)={}".format(a))
    plt.ylabel("blocksize of repetition coding")
    plt.xlabel("Probability of bit mismatch")
    plt.yticks([1.5,2.5,3.5,4.5], ["(no rep) 1","2","3","4"])
    plt.show()
        
