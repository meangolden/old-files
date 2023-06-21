# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 19:01:24 2021
Description: lighweight RFID protocol V2 SECTION

@author: cp17593
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.interpolate import make_interp_spline, BSpline
from math import pi, factorial
from scipy.signal import argrelextrema
def spatCorr(model,D,l=1):
    '''Inputs:  model = 
                D= distance vector in m  
                l= wavelength, for the default value l=1
                    the distance is normalised to the wavelength
        Output: spatial correlation in the 3D diffuse field with Rayleigh
                channels against D (Type:list). 
                '''
    assert model in ["2D","3D"]
    if model =="2D":
        corr = [jv(0,2*pi*d/l) for d in D ]
    elif model =="3D":
        corr = [np.sin(2*pi*d)/(2*pi*d) if d>0 else 1 for d in D]
    else:
        print("check model type")
    return corr
    
def nCr(n,r):
    return factorial(n)/(factorial(r)*factorial(n-r))

def probFraud(M,epsilon):
    ''' Returns the probability of solo distance fraud when the size
        of the link signature is M and the decision threshold is epsilon'''
    assert isinstance(M, int)
    summ=0
    epsilon= (epsilon+1)/2
    for k in range(int(M*epsilon)+1):
        summ = summ + (-1)**k*nCr(M,k)*(M*epsilon-k)**M 
    return 1-(1/factorial(M))*summ
    

if __name__ == "__main__":
    #choice = input("Press 2 for 2D \n Press 3 for 3D.\n")
    #model = "2D" if choice=="2" else "3D" 
    model = "3D"
    
    D = np.linspace(0,6,300)   
    M = range(1,30)
    
    corr = spatCorr(model,D,l=1)
    indx = argrelextrema(np.array(corr),np.greater)[0][0]
    min_threshold = corr[indx]
    max_distance = D[indx]
    
    print("maxdistance= ",max_distance )
    print("min_threshold= ",min_threshold)
    
    #E = np.linspace(min_threshold,0.8,200)
    E = np.linspace(0,0.8,200)
    P =[probFraud(m,e) for e in E for m in M]
    Pm = np.reshape(P,(len(E),len(M)))
    Em, Mm = np.meshgrid(M,E)

    
    fig = plt.figure(figsize=(4,4))
    ax = Axes3D(plt.gcf())
    my_cmap = plt.get_cmap('coolwarm') 
    surf = ax.plot_surface(Em, Mm, Pm, cmap=my_cmap)
    ax.set_ylabel("$\epsilon$",fontsize=12)
    ax.set_xlabel("$M$",fontsize=12)
    ax.set_zlabel("$P$(missed detection)",fontsize=12)
    fig.colorbar(surf, shrink=0.6, aspect=10)

    plt.show()


