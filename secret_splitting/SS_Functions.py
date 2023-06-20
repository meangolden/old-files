#!/usr/bin/env python3

# Chrys 2020-04-02

# Standard library
from math import *

# PyPI (Python Package Index) library
import numpy as np 
from numpy import linalg as LA
import scipy.linalg as sla

def isIntegral(x):
    ret = (isinstance(x, int)) or \
          (np.issubdtype(x, np.integer))
    return ret

def isFloating(x):
    ret = (isinstance(x, float)) or \
          (np.issubdtype(x, np.floating))
    return ret

def isReal(x): #I do not like this function. It doesn't tell you if a number is complex or real! Change or discard.
    ret = isIntegral(x) or \
          isFloating(x)     
          
def isNumber(x):
    ret = (isinstance(x, (float, int, complex))) or \
          (np.issubdtype(x, np.number))
    return ret

def assertIntegral(x):
    assert isIntegral(x), (type(x), x)

def assertNumber(x):
    assert isNumber(x), (type(x), x)

def channelMatrix(N, NA, sigma_2):
    '''Generates a matrix of size NAxN whose entries follow the CN(0,sigma_2) distribution.
    Reminder: Τhe inputs for the function normrdn(μ,σ) are the mean and
    the standard deviasion. NOT the variance.
    '''
    assertIntegral(NA)
    assertIntegral(N)

    sigma = sqrt(sigma_2/2)
    H = np.random.normal(0,sigma,(N,NA)) + np.random.normal(0,sigma,(N,NA)) * 1j
    return H

def noeveChannCap(P, NA):
    Cb_ = 0.0;
    for i in range(0, 1000):
        H = np.random.normal(0,sqrt(.5),(1, NA)) + \
            np.random.normal(0,sqrt(.5),(1,NA)) * 1j
        Cb_ = Cb_ + max(0, P*np.log2(1 + np.linalg.norm(H)**2))
        #l = max(np.linalg.eigvalsh(H @ np.conjugate(H.T))) #same result as line above
        #assert isNumber(l)
        #Cb_ = Cb_ + np.log2(1+ l)
    return Cb_/100

def channelStats(xe, ye, alpha = 4):
    '''Inputs: (xe, ye): are the cartesian coordinates of the eavesdropper's location
               when the base stations are placed at A1(1, 0) and A2(-1, 0).
               alpha: the pathloss exponent
       Returns: variances noisePowere1 and noisePowere2 of the channel coefficients ~ CN(0,noisePowerei)'''
    #A1, A2 = (1, 0), (-1, 0)
    dist1, dist2 = np.sqrt((xe-1)**2+(ye)**2), np.sqrt((xe+1)**2+(ye)**2)

    return 1/dist1**alpha, 1/dist2**alpha

def ctranspose(H):
    '''Returns the conjugate transpose of matrix H'''
    return np.conjugate(H.T)

def gram(H):
    '''The gram of a matrix H is equal to the product of itself and its Hermitian
    return. Reminder: Bob's and Eve's matrices are of shape (Nb, Na), (Ne, Na).
    For the other way around, change the order you multiply.'''
    
    return(ctranspose(H) @ H)
   
def isSemidefinite(H): #haven't been used but I may need it. 
    '''Input: a square matrix. Returns: True if H is semidifinite, False, otherwise.'''
    assert np.shape(H)[0] == np.shape(H)[1], print("Input matrix must be square.")
    E= np.linalg.eigvalsh(H)
    return all( e > 0 for e in E)
  
def isNonDefinite(H):
    '''Returns True or False whether matrix H is non definite or not. 
       Recall that a matrix is non definite if there exists at least one
       positive eigen value and at least one negative value.'''
    E= np.linalg.eigvalsh(H)
    return any( e > 0 for e in E) and any(e < 0 for e in E)

def expmimomeSplitCs(sigma1, sigma2, P, NA1, NA2, NB, NE, nIter):
    assert NA1 == 2, ("channel model must me 2-1-1, or 2-2-2, or 2-1-2")
    assert NA2 == 2, ("channel model must me 2-1-1, or 2-2-2, or 2-1-2")
    assert NB < 3, ("channel model must me 2-1-1, or 2-2-2, or 2-1-2")
    assert NE < 3, ("channel model must me 2-1-1, or 2-2-2, or 2-1-2")
    Cs_ = 0.0 #initialisation
    for i in range(nIter):
        Hb1 = channelMatrix(1, NA1, sigma1) #Bob is a single antenna node
        Hb2 = channelMatrix(1, NA2, sigma1)
        He1 = channelMatrix(NE, NA1, sigma2)
        He2 = channelMatrix(NE, NA2, sigma2)
        c1 = mimomeCs(Hb1, He1, P)
        c2 = mimomeCs(Hb2, He2, P)
        Cs_ = Cs_ + splitCs(c1, c2)
        #print(c1, c2, splitCs(c1, c2))
        assert Cs_ >= 0.0
    return Cs_ / nIter

def eveSnrMaxima(eveSNRs):
    ''' This function takes as input a set of eavesdroppers' SNRs. 
    eveSNRS is a K x 2 matrix; The kth row contains the SNRs between the kth
    eavesdropper A1, and the SNR of the same eavesdropper with A2.
    the function lodivide them into two subsets. The first subset contains the eavesdroppers
    that have a better channel to A2. It returns the maximum SNR from each 
    subset.'''
    assert np.shape(eveSNRs)[1] == 2
       
    subset1 = [] # gathers the smallest SNRs for the channel of A1
    subset2 = []
    l = len(eveSNRs)
    for i in range(l):
        if eveSNRs[i][0] < eveSNRs[i][1]:
            subset1.append(eveSNRs[i][0])
        else:
            subset2.append(eveSNRs[i][1])
    
    if subset1 == []:
        assert len(subset2) == l
        maxSNR1 = 0
        maxSNR2 = np.max(subset2)
    elif subset2 == []:
        maxSNR2 = 0
        maxSNR1 = np.max(subset1)
    else:
        maxSNR1, maxSNR2 = np.max(subset1), np.max(subset2)
    
    return maxSNR1, maxSNR2
    
# Functions that calculate Secrecy Capacity

def splitCs(C1, C2): 
    ''' MUST REVIEW'''    
    
    return max(C1, C2)
  
def NaTwoCs(Hb, He, P, noisePower=1): #seems to work well
    ''' Returns the secrecy capacity for the case when Alice has two antennas, 
        gram(Hb) > gram(He), Bob and Eve have the same power noise, and not more
        than 2 antennas each. '''

    assert Hb.shape[1] == 2 and He.shape[1] == 2, "hi" #Alice has two antennas.
    assert Hb.shape[0] < 3 and He.shape[0] < 3, "hello" #The receivers have up to 2 antennas.
    assert isNonDefinite(gram(Hb) - gram(He)) == True, 1
    
    assert gram(He)[1][0] == np.conjugate(gram(He)[0][1])
    assert gram(Hb)[1][0] == np.conjugate(gram(Hb)[0][1])
    
    
    rho = P / noisePower
    a1, b1, c1 = gram(He)[0][0], gram(He)[0][1], gram(He)[1][1]
    a2, b2, c2 = gram(Hb)[0][0], gram(Hb)[0][1], gram(Hb)[1][1]
    
    assert all(np.isscalar(x) for x in (a1, a2, b1, b2, c1, c2))
    assert all(np.iscomplexobj(x) for x in (a1, a2, b1, b2, c1, c2))
    
    
    A = a1 * c2 + a2 * c1 - b1 * np.conjugate(b2) - np.conjugate(b1) * b2 + \
        (a1 + a2 + c1 + c2) / rho + 2 / rho**2
    #print("A", A)
    prod1 = (a2 * c1 - b1 * np.conjugate(b2) + (a2 + c1) / rho + 1/rho**2) 
    prod2 = (a1 * c2 - np.transpose(b1) * b2 + (a1 + c2) / rho + 1/rho**2)
    prod3 = (b2 * c1 - b1 * c2 + (b2 - b1) / rho)
    prod4 = (a1 * np.conjugate(b2) - a2 * np.conjugate(b1) + (np.conjugate(b2) - np.conjugate(b1)) / rho)
    B1 = prod1 * prod2
    B2 = prod3 * prod4
    #print("B1 = ", B1)
    #print("B2 = ", B2)
    D = A**2 - 4 * (B1 - B2)
    #print("D = ", D)
    numerator = A + np.sqrt(D)
    denominator =  2 * ((1/rho + a1)*(1/rho + c1) - LA.norm(b1)**2)
    Cs = np.log2(numerator / denominator)
    #assert round(np.imag(Cs), 5) == 0
    
    return np.real(Cs)

def misoseCs(Hb, He, P, noisePower=1): #this formula gave similar results with other methods but not the same
    '''Inputs: Bob's channel vector, Eve's channel vector, transmit power, noise
               power. It is assumed that both receivers have the same noise power
       Returns: The secrecy capacity for the MISOSE case. '''
    Hb = np.array(Hb)
    He = np.array(He)
    assert Hb.shape[0] == 1 #Bob is a single antenna node 
    assert He.shape[0] == 1 #Eve is a single antenna node
    rho = P / noisePower
    a = 1 + rho * LA.norm(He)**2
    b1 = 2 + rho * LA.norm(He)**2 + rho * LA.norm(Hb)**2
    b2 = LA.norm(Hb)**2 * LA.norm(He)**2
    b3 = (LA.norm(ctranspose(Hb.T) @ He.T ))**2
    b = b1 + rho**2 * (b2 - b3)
    
    c = 1 + rho * LA.norm(Hb)**2
    Cs = np.log2((b + np.sqrt(b*b - 4*a*c))/(2*a))
    return np.log2((b + np.sqrt(b*b - 4*a*c))/(2*a))
    
def misomeCs(Hb, He, P): # H: Nr x Nt
    assert Hb.shape[0] == 1 #Bob is a single antenna node
    W1 = ctranspose(Hb) @ Hb
    W2 = ctranspose(He) @ He
   
    assert W1.shape == W2.shape

    I = np.identity(np.size(W1,0))
    A, B = I + P*W1, I + P*W2
    ctranspose(B)
    #E= np.linalg.eigvalsh(np.linalg.inv(B) @ A)
    E, _= sla.eig(A, B) #perfomrs the same as above line. Probably more stable. 
    e = np.amax(E)
    #print("max eigen vaue:", e)
    assert np.isscalar(e)
    Cs = max(np.log2(e), 0.0)
    #print("Ctrad = ", Cs)
    #assert np.imag(Cs) < 0.1**2
    return np.real(Cs)


#Expectation
def expmisomeSplitCs(sigma1, sigma2, P, NA1, NA2, NE, nIter):

    Cs_ = 0.0 #initialisation
    for i in range(nIter):
        Hb1 = channelMatrix(1, NA1, sigma1) #Bob is a single antenna node
        Hb2 = channelMatrix(1, NA2, sigma1)
        He1 = channelMatrix(NE, NA1, sigma2)
        He2 = channelMatrix(NE, NA2, sigma2)
        if LA.norm(He1) > LA.norm(He2):
            c1 = log2(1 + LA.norm(Hb1))
            c2 = misomeCs(Hb2, He2, P)
        else:
            c2 = log2(1 + LA.norm(Hb2))
            c1 = misomeCs(Hb1, He1, P)
            
        Cs_ = Cs_ + splitCs(c1, c2)
        #print(c1, c2, splitCs(c1, c2))
        assert Cs_ >= 0.0
    return Cs_ / nIter

def expNaTwoCs(sigmaE_2, P, NA, NE, nIter):

    Cs_ = 0.0 #initialisation
    for i in range(nIter):
        Hb = channelMatrix(1, NA, 1) #Bob is a single antenna node
        He = channelMatrix(NE, NA, sigmaE_2)
        c = NaTwoCs(Hb, He, P)
        #cb = np.log2(1 + P*np.linalg.norm(Hb)**2) #when there is no eavesdropper MRT maximises Bob's SNR
        Cs_ = Cs_ + max(0 , c)

    assert isinstance(Cs_, float), (type(Cs_), Cs_)
    return Cs_ / nIter

def expmisoseCs(sigmaE_2, P, NA, NE, nIter):

    Cs_ = 0.0 #initialisation
    for i in range(nIter):
        Hb = channelMatrix(1, NA, 1) #Bob is a single antenna node
        He = channelMatrix(NE, NA, sigmaE_2)
        c = misoseCs(Hb, He, P)
        #cb = np.log2(1 + P*np.linalg.norm(Hb)**2) #when there is no eavesdropper MRT maximises Bob's SNR
        Cs_ = Cs_ + max(0 , c)

    assert isinstance(Cs_, float), (type(Cs_), Cs_)
    return Cs_ / nIter

def expmisomeCs(sigmaE_2, P, NA, NE, nIter):

    Cs_ = 0.0 #initialisation
    for i in range(nIter):
        Hb = channelMatrix(1, NA, 1) #Bob is a single antenna node
        He = channelMatrix(NE, NA, sigmaE_2)
        c = misomeCs(Hb, He, P)
        #cb = np.log2(1 + P*np.linalg.norm(Hb)**2) #when there is no eavesdropper MRT maximises Bob's SNR
        Cs_ = Cs_ + max(0 , c)

    assert isinstance(Cs_, float), (type(Cs_), Cs_)
    return Cs_ / nIter
  
#Artificial Noise Generation

def transmittingwithAN(Hb):
    '''Input: Bob's channel.
       Outputs: principalEigValue, beamformer (same as MRT beamformer because 
       of MISOME case), noisevector, ANcovariance.
       
       When there is no knowledge of the eavesdropper's channel, AN generation is often 
       employed in order to decrease her channel.
       For theory advise paper: PL secrecy of MIMO communication in the presence
       of a Poisson random filed of eavesdroppers.'''
    assert np.shape(Hb)[0] == 1, "Bob must be a single antenna (MISOME)"
    
    NA = np.shape(Hb)[1]
    eigValues, eigVectors = LA.eigh(gram(Hb))
    eigValues = np.round(eigValues.real, 5) + np.round(eigValues.imag, 2) * 1j
    eigVectors = np.round(eigVectors.real, 5) + np.round(eigVectors.imag, 2) * 1j
        
    principalEigValue = eigValues[-1]
    assert np.sum(eigValues) == principalEigValue # we're expecting a single
    # positive eigenvalue. That's because the rank of H (and rank(H^HH)) is 1.
    beamformer = eigVectors[-1]
    
    summ = np.sum(eigVectors, axis=0) - beamformer
    assert len(summ) == NA
    noisevector = sqrt(1/(NA - 1)) * summ
    assert NA > 1 , "this method does not work for NA = 1"
    ANcovariance = gram(np.asmatrix(noisevector))
    ANcovariance = ANcovariance.real
    assert round(np.matrix.trace(ANcovariance).item(0,0), 1) == 1, round(np.matrix.trace(ANcovariance).item(0,0), 2)
    assert np.shape(ANcovariance) == (NA,NA)
    
    return np.real(principalEigValue), beamformer, noisevector, ANcovariance
    
def SNRBob(phi, P , maxEigenValue, noisePower = 1):
    '''When there is no knowledge of the eavesdropper's channel, AN generation is often 
       employed in order to decrease her channel. Alice transmits 
       sqrt(phi)Ptd +sqrt(1-phi)Pn. t is the beamformer: the eigen value that 
       corresponds to the largest eigen value of Hb^HHb. In this case, i.e., in this
       MISOME case, t1 is simply the MRT beamformer. 
     
       For theory advise paper: PL secrecy of MIMO communication in the presence
       of a Poisson random filed of eavesdroppers. '''
    assert phi <= 1
    infoPower = phi * P # the trasnsmit power allocated fot the information signal. 
    snr = (infoPower * maxEigenValue) / noisePower
    assert snr.imag < 0.1 ** 4
    return snr.real
    
def SNREve(phi, P, beamformer, ANcovariance, He, noisePower=1): 
    '''assertion error: some time it gives a complex SNR (the imag. of 
       which is not insignificant.)
       When there is no knowledge of the eavesdropper's channel, AN generation is often 
       employed in order to decrease her channel. Alice transmits 
       sqrt(phi*P)td +sqrt((1-phi)*P)n. n is the noise vector which is orthogonal to 
       the bemaformer used for the the data signal. Eve is assumed to know the beamformer, 
       the covariance of noise vector, her channel, and phi. She employes receive
       breamforming to increase her SNR (worst case scenario).
       For theory advise paper: PL secrecy of MIMO communication in the presence
       of a Poisson random filed of eavesdroppers. '''
    assert phi <= 1
    t1 = beamformer
    Cn = ANcovariance
    NE = np.shape(He)[0]
    IntNoise = (phi) * P * He @ Cn @ ctranspose(He) + noisePower * np.eye(NE)
    SNRe = (1 - phi) * P * ctranspose(t1) @ ctranspose(He) @ LA.inv(IntNoise)@ He @ t1
    SNRe = SNRe.real
    #assert SNRe.imag < 0.1 ** 2
    assert SNRe.shape == (1, 1)
    return SNRe.item(0,0)
    
def SNREve1(phi, P, beamformer, noisevector, He, noisePower=1):
    ''' According to paper PLS with AN: Secrecy Capacity and optimal 
        Tx Power allocation.    '''
    assert He.shape[0] == 1, "it has been seen that this works only when Eve \
                                is a single antenna node"
    numerator = phi * P * LA.norm(He @ beamformer)**2
    denom = (1 - phi) * P * LA.norm(He @ noisevector)**2 + noisePower
    snr = numerator/denom
    #assert snr.shape == (1,1) 
    return np.reshape(snr, -1)[0]
    
    
def SNRBob1(phi, P , Hb, noisePower=1):
    '''It perrfoms the same as SNRBob. Double checked!
       When there is no knowledge of the eavesdropper's channel, AN generation is often 
       employed in order to decrease her channel. Alice transmits 
       sqrt(phi)Ptd +sqrt(1-phi)Pn. t is the MRT beamformer. '''
    assert phi <= 1
    infoPower = phi * P # the trasnsmit power allocated fot the information signal. 
    snr = (infoPower * LA.norm(Hb)**2) / noisePower
    assert snr.imag < 0.1 ** 4
    return snr.real
    
    
def RsANbmfEve(phi, p, Hb, He):
    maxEig, beamformer, noisevector, ANcovariance = transmittingwithAN(Hb)
    SNRb = SNRBob(phi, p , maxEig)
    SNRe = SNREve(phi, p, beamformer, ANcovariance, He)
    RsANwithBMFatEve = max(0, (log2(1 + SNRb) - log2(1 + SNRe)))
    return RsANwithBMFatEve
    #breakpoint()
    
def RsANbmNe1(phi, p, Hb, He):
    assert He.shape[0] == 1, "this method works only for the case when NE = 1"
    maxEig, beamformer, noisevector, ANcovariance = transmittingwithAN(Hb)
    SNRb = SNRBob(phi, p , maxEig)
    SNRe1 = SNREve1(phi, p, beamformer, noisevector, He)
    RsAN = max(0, (log2(1 + SNRb) - log2(1 + SNRe1)))
    return RsAN
    
#def RsANnobmfEve(phi, p, Hb, He):
#    maxEig, beamformer, noisevector, ANcovariance = transmittingwithAN(Hb)
#    SNRb = SNRBob(phi, p , maxEig)
#    SNREve = 
    
    