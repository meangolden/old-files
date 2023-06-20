# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 15:34:17 2020

@author: chrys
"""

#Standard library
from math import *
import matplotlib.pyplot as plt
import numpy as np

def nCr(n,r):
    numer = 1
    for i in range(r):
        numer = numer * (n - i)
    assert numer % factorial(r) == 0
    return int(numer / factorial(r))
        
    return c

def pr_rej_Alice(a,l,t):
    '''Inputs: a: Probability of a single bit mismatch
               l: Lenght of sequence
               t: Threshold: Bob accepts Alice if the hamming distance
                  between his seq. and Alice's seq is not greater to t.
       Output: The probability of rejecting Alice. '''
    assert a <= 1
    assert t <= l
    assert isinstance(l,int)
    assert isinstance(t,int)
    s = 0.0
    for m in range(0,t+1):
        s += nCr(l, m) * a**m * (1-a)**(l-m)
        
    pr = round(1-s,5)
    assert(pr >= 0)
    return pr

def pr_acc_Eve(l,t):
    #CHECK IF THIS IS CORRECT
    '''Inputs: l: Lenght of sequence
               t: Threshold: Bob accepts Alice (hence rejects Eve)
                  if the hamming distance between his seq. and her seq is not
                  greater to t.
       Output: The probability of accepting Eve. '''
    s = 0
    for fg in range(0,t+1):  #false guesses
        s += nCr(l, fg) * 0.5**l
    assert round(s,5) <=1 , s
    return s

def probSameBit(eb,ea):
    '''Bob and Alice observe Randy's binary trasmissions. The channel between
    Bob and Randy is a binary symmetrical channel with cross over prob eb.
    Similarly, we define ea for Alice's channel. 
    Output: the probability that Bob and Alice observe the same bit (in one
    channel use.)'''
    p = eb*ea + (1-eb)*(1-ea)
    assert p <= 1
    assert p >=0
    return p

def probDifBit(eb,ea):
    '''Bob and Alice observe Randy's binary trasmissions. The channel between
        Bob and Randy is a binary symmetrical channel with cross over prob eb.
        Similarly, we define ea for Alice's channel. 
        Output: the probability that Bob and Alice observe different bit (in one
        channel use.)'''
        
    return 1 - eb*ea - (1-eb)*(1-ea)


def probSameKeybit(m,e):
    '''Inputs: m: The blocksize of repetition coding;
               e: The prob of dif Bit.
       Output: the probability that repetition coding has been successful for 
               the block, resutling in an identical key bit for Bob and Alice.
    '''
    p = 0
    for i in range(int(m/2)+1, m+1):
        p = p + nCr(m,i) * (1-e)**i * (e)**(m-i)
        assert p <= 1
        assert p >= 0
    return p
    
def probdifKeybit(m,e):
    '''Inputs: m: The blocksize of repetition coding;
               e: The prob of differnt Bit.
       Output: the probability that repetition coding has been unsuccessful for 
               the block, resutling in an identical key bit for Bob and Alice.
    '''
    p = 0
    for i in range(int(m/2)+1, m+1):
        p = p + nCr(m,i) * (e)**i * (1-e)**(m-i)
    assert p <= 1
    assert p >= 0
    return p

def cantDecide(m,e):
    c = round(1 - probSameKeybit(m,e) - probdifKeybit(m,e),3)
    assert c >= 0, "prob of mismatch is assumed to be less than .5"
    
    return c
    
    

 