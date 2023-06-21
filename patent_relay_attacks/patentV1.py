# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:42:12 2021

@author: cp17593
"""
# mport numpy as np
from numpy import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from math import *

N = range(2,32) # n in N is the size of the channel sequence
K = range(1,34) # k in K is the size of the output of the hash function
L = range(1,34) # l in L is the size of the nonce sent by the the verifier
E = 0.3 #percentage threshold
corr_coeff_attack = []

#A relay attack is successful if the verifier and the the prover, located far apart, 
# experience the similar channel, i.e., they attain similar channels. The hamming di

def nCr(n,r):
    return factorial(n)/(factorial(r)*factorial(n-r))


Pr_success_relay = []
i = 0
for n in N:
    e = int(E*n) #min(N)#int(E*n) # e is the threshold, i.e. allowed mismatches in the channel sequence
    summ = 0
    for k in range(e+1):
        summ = summ + nCr(n,k)*0.5**n
    print(summ)   
    print(i)
    Pr_success_relay.append(summ)
    i  = i + 1
    
# In a replay attack, the adversary observes a common channel.I.e., she knowns the channel sequence.
#But without the secret key, she doesn't know the hash function (of size K).
    

average_no_guesses = [2**(k-1) for k in K] # simply guessing.
average_no_queries = [2**(l-1) for l in L] #brute force: the attacker asks for "link signatures"

#fig, ax = plt.subplots()
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels)
#plt.plot(Pr_success_relay, linewidth = "4", label = "e = 12 jhk")
plt.legend()
fig, ax = plt.subplots()
#plt.plot(K, average_no_guesses)
plt.show()
plt.plot(average_no_queries,linewidth = "3")
#plt.yscale('log')
