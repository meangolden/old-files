# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 21:34:43 2020

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
    #it gives the probability of generating an identical key after repetition coding.
  
    
    #variables
    lngth = 300 #length of correlated sequences
    M = [1, 2, 3, 4, 5, 6,7,8] #blocksize of repetition coding (M=1 is the trivial case. No repetition.
    L = [int(lngth/m) for m in M] # the nunber of blocks (key size may be smaller)
    
    eb = 0.1#cross over probability for Bob
    Ea = np.linspace(0.1,0.3,num=3) #cross over probability for Alice (many cases)
    
    Ed = [myf.probDifBit(eb,ea) for ea in Ea]
    Es = [myf.probSameBit(eb,ea) for ea in Ea]
    
    #probability of a block leading to a different key bit. 
    Pd = [[round(myf.probdifKeybit(m,e),2) for e in Ed] \
          for m in M]
    #probability of a block leading to the same key bit
    Ps = [[round(myf.probSameKeybit(m,e),2) for e in Ed] \
          for m in M]
    noBlk = [int(lngth/m) for m in M]

    #probability that there exists a "bad" block leading to diff key bit(s), 
       #thus "wrong" key.    
    def prWrongKey(m, e):
        noBlck = int(lngth/m)
        return 1-(1-myf.probdifKeybit(m,e))**noBlck
    
    Pwk = [[round(prWrongKey(m, e), 2) for e in Ed] \
           for m,b in zip(M,noBlk)]
        
    ExpMisMatchs = [[round(myf.probdifKeybit(m,e)*nb,3) for e in Ed] \
                     for (m,nb) in zip(M, noBlk)]
        
    ExpLength = [[round(nb - myf.cantDecide(m, e)*nb,3) \
                  for e in Ed] for (m,nb) in zip(M, noBlk)]
    #graph
    #plt.plot(Ea,A)
    #plt.plot(Ea,B)
    #plt.plot(Ea,C)
    #plt.xlabel("Prob of channel estimation error for Alice")
    #plt.title("Probability of generating a perfect matching key. Prob of channel estimation error for Bob is set to 0.1")
    #plt.yticks([1.5,2.5,3.5,4.5], ["(no rep) 1","2","3","4"])
    #plt.show()
        
