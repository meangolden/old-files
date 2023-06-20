# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:34:15 2020

@author: chrys
"""

#Standard library
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import signal

def manchester(binary_data):
    for i in binary_data:
        assert (i == 0 or i == 1), "this is not binary"
    
    m = []
    for b in binary_data:
        if b == 0:
            m.append(1)
            m.append(0)
       
        else:
            m.append(0)
            m.append(1)
             
    plt.plot(m, marker='', color='blue', drawstyle='steps-post')
    for i in range(len(binary_data)):
        plt.text(2*i+1, 1.1 , str(binary_data[i]), fontsize=12)
        
    return(m)

def miller(binary_data):
    if binary_data[0] == 0:
        m = [0,1]
    else:
        m = [1,0]
    binary_data.pop(0)
        
    for b in binary_data:
        if b == 0:
            m.append(m[-1])
            m.append(m[-1])
       
        else:
            m.append(m[-1])
            m.append((m[-1]+1) % 2)
            
    plt.plot(m, marker='', color='red', drawstyle='steps-post')
    for i in range(len(binary_data)):
        plt.text(2*i+1, 1.1 , str(binary_data[i]), fontsize=12)
            
    return m
    
    
    

