# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 13:51:32 2021

@author: chrys
"""

#Plots the probability of missed detection for a fixed upper bound of 
#flase alarm probability, against the probability of mismatch. Recall that 
#Pro of missed detection does not change directly with p, but it changes with t,
#and t is function of p when there is a threshold on the false alarm probability.
#standard libary
from math import *
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
#mylibrary
import myfunc as myf
#W R O N G
# keep only for reference to graphs
if __name__ == "__main__":
    
    
    #fixed variables
    P = np.linspace(0,0.35,20) #probability of bit mismatch for Bob and Alice
    l = 200 #legth of correlated sequences
    
    
    i = 0
    missedDet = 8*np.ones(len(P))
    falseAlarm = 8*np.ones(len(P))
    thr = l * np.ones(len(P))
    tminimum = l * np.ones(len(P))
    fa = 0.05 #threshold for the false alarm probability
    for p in P:
        while thr[i]>=0:
            prFA = myf.pr_rej_Alice(p,l,int(thr[i]))
            if prFA <= fa:
                falseAlarm[i] = prFA
                missedDet[i] = myf.pr_acc_Eve(l,int(thr[i]))  #value will be ovewritten. Last value
                #corresponds to the min threshold for which P(reg Alice)<=0.01
                tminimum[i] = int(thr[i])
            thr[i] -=1
        i = i + 1
        
    #print("missedDet",missedDet)
    #print("falseAlarm",falseAlarm)
    #print("tminimum",tminimum)

 
    #graphs
    #plt.plot(P,falseAlarm, linewidth=3, markersize=18)
    plt.plot(P, missedDet, linewidth=2, markersize=18)
   # plt.legend(["rejecting Alice (false alarm)","accepting Eve (missed detection)"])
    plt.title("Pr(missed detection) when Pr(false alarm)<={} and length={}".format(fa,l))
    plt.ylabel("p.d.f.")
    plt.xlabel("q")
    plt.grid(color='black', linestyle='-', linewidth=1)
    #plt.show()
    for i in range(len(P)):
        if i % 6 == 0:
            plt.text(P[i], missedDet[i], r't={}'.format(int(tminimum[i])),fontsize=14)
            plt.scatter(P[i], missedDet[i], c=(0,0,0))
    



