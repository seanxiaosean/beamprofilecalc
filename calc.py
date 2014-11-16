# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 11:11:01 2014

@author: xiao

Basic calculation for 1030nm laser system.

check the unit in propagation calc.

default unit is mm.

all beam diameter/radius is in d63 standard
"""
import numpy as np
import matplotlib.pyplot as plt
  
pi = 3.141592653589793

# system wavelength
wl = 1.03e-6
# default system M^2 value
M2 = 1.2
# default propagation distance between the telescope and the coupling lens!
pd = 4000


def bpp(M2):
    '''
    unit: beam radius:mm angle:rad
    bpp(1) : return diffraction limit lambda/pi   
    '''        
    BPP = M2 * (wl*1.e3/pi)

    return BPP
    
def gaussianlens(z,w,f):
    '''
    all units are in mm.    
    
    '''
    w1 = np.sqrt(w**2/((1-z/f)**2+(pi*w**2/(wl*1e3)*f)**2))
    
    z1 = f - (f-z)/((1-z/f)**2+(pi*w**2/(wl*1e3)*f)**2)
    
    return z1,w1
     
# define a function to calculate the beam diameter for a random z
   
    
    
def gaussiantelescope(z,w,f1,f2,d):
          
    z1prime,w1 = gaussianlens(z,w,f1)
     
    z2 = f1+f2+d-z1prime
     
    z2prime,w2 = gaussianlens(z2,w1,f2)

    theta2 = bpp(M2)/w2
    
    w_after_chirpmirrors = gaussianbeam(w2,pd-z2prime)
    
    return z2prime,w2,theta2,w_after_chirpmirrors
    
    
     

def gaussian(w,z):
    
    wz = w*np.sqrt(1+((wl*1e3)*z/(pi*w**2))**2)
    
    return wz

def gaussianbeam(w,z):
    '''
    for given pulse parameters, calculate the beam radius at certain place
    '''
    w2 = w*np.sqrt(1.+((wl*1e-3)*z/(w**2*pi))**2)

#    r = z*(1.+((wl*1e3)*z/(w**2*pi))**-2)**2    
    
    return w2

def plotgaussian(w):
    '''
    Plot gaussian beam profile within +-3z_r by given wrist radius
    '''
    z_r = pi*w**2/1.03e-3
    
    z = np.linspace(-3*z_r,3*z_r,100)
    w2 = w*np.sqrt(1.+((wl*1e3)*z/(w**2*pi))**2)

    plt.plot(z,w2,color='blue')
    plt.hold(True)
    plt.grid(True)
    plt.plot(z,-w2,color='blue')
    plt.axvline(x=z_r,color='black',linestyle="--")
    plt.axvline(x=-z_r,color='black',linestyle="--")
    plt.hold(False)
   
def bppfiber(NA,d):  
    '''
    NA and MFD!
    d unit in um
    BPP unit: mm rad
    '''
    BPP = natoangel(NA) * d/2 * 1e-3
    return BPP
   
   
def natoangel(NA):

    theta = np.arcsin(NA)     
    
    return theta
    
def angel(w):
    '''
    unit:mrad
    '''
    theta = (wl*1e3)/(pi*w)*1000

    print str("%.4f" %theta) + ' mrad'
    return theta    
    
    

def coupling(w,MFD):
    
    f = lensfocus(2*w,MFD)
        
    z,w = gaussianlens(f,w,f)
    
    return z,w,f


def lensfocus(Dbeam63,MFD):
    '''
    Dbeam is the diameter of D63(mm)
    f in mm
    MFD is the Mode Field Diameter of fiber(micron)
    '''
    f = Dbeam63*np.sqrt(2)*MFD*0.7625
    
    return f
    
