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
#W R O N G
# keep only for reference to graphs
if __name__ == "__main__":
    
    
    #fixed variables
    eb = 0.1 #flip probability for Bob
    ea = 0.2 #flip probability for Alice
    l = 100 #legth of correlated sequences
    
    a = eb*(1-ea) + (1-eb)*ea 
    
    falseAlarm = [myf.pr_rej_Alice(a,l,t) for t in range(int(l*.75))]
    missedDet = [myf.pr_acc_Eve(l,t) for t in range(int(l*.75))]
 
    #graphs
    plt.plot(range(int(l*.75)),falseAlarm, linewidth=3, markersize=18)
    plt.plot(range(int(l*.75)),missedDet, linewidth=3, markersize=18)
    plt.legend(["rejecting Alice (false alarm)","accepting Eve (missed detection)"])
    plt.title("Pr(bit mismatch)={}".format(a))
    plt.ylabel("p.d.f.")
    plt.xlabel("t (#allowed mismatches)")
    plt.grid(color='black', linestyle='-', linewidth=1)
    plt.show()
        
