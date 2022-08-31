import numpy as np
import os
import sys
import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

#################

def lisa_noise(ff):

    Larm = 2.5e9 # m
    fstar = 19.09e-3 # mHz

    pn = lisa_sn(ff, Larm, fstar)

    sc = lisa_sc(ff)

    rn = 3.0 / 10.0 / (1.0 + 0.6 * (ff / fstar) ** 2.0)

    yvalue = pn / rn + sc

    return yvalue

def lisa_sn(ff, Larm, fstar):

    P_oms = (1.5e-11) ** 2 * (1. + (2.0e-3 / ff) ** 4)

    P_acc = (3.0e-15) ** 2 * (1. + (0.4e-3 / ff) ** 2) * (1. + (ff / (8.0e-3)) ** 4)

    Pn = (P_oms + 2.0 * (1.0 + np.cos(ff / fstar) ** 2) * P_acc / (2.0 * np.pi * ff) ** 4) / Larm ** 2
    
    return Pn

def lisa_sc(ff):

    AA = 9e-45
    alpha = 0.138
    beta = -221.0
    gamma= 1680
    kappa = 521
    fk = 0.00113

    sc = 1. + np.tanh(gamma*(fk - ff))
    sc *= np.exp(-ff ** alpha + beta * ff * np.sin(kappa * ff))
    sc *= AA * ff ** (-7.0 / 3.0)

    return sc

def ligo_noise(ff):

    x = ff / 245.4

    yvalue = 1e-48 * (0.0152 / x ** 4. + 0.2935 * x ** 2.25 + 2.7951 * x ** 1.5 - 6.5080 * x ** 0.75 + 17.7622)

    return yvalue

def decigo_noise(ff):

    x = ff / 7.36

    yvalue = 6.53e-49 * (1. + x ** 2.) + 4.45e-51 / ff ** 4. / (1. + x ** 2.) + 4.94e-52 / ff ** 4.

    return yvalue

def et_noise(ff):

    data = np.loadtxt('et_noise.txt')

    x = data[:,0]
    y = data[:,1]

    (Idx,) = np.where(x > ff)
    i = Idx[0] - 1
    yvalue = (y[i+1] - y[i]) / (x[i+1] - x[i]) * (ff - x[i]) + y[i]

    yvalue *= yvalue

    return yvalue

def hc(ff, m1, m2, z):

    mchirp = (m1 * m2) ** (3.0 / 5.0) / (m1 + m2) ** (1.0 / 5.0)
    eta = m1 * m2 / (m1 +  m2) ** 2.

    a0, b0, c0 = 2.9740e-1, 4.4810e-2, 9.5560e-2
    f0 = (a0 * eta ** 2 + b0 * eta + c0) / (m1 + m2) * 2e5 / np.pi 
    a1, b1, c1 = 5.9411e-1, 8.9794e-2, 1.9111e-1
    f1 = (a1 * eta ** 2 + b1 * eta + c1) / (m1 + m2) * 2e5 / np.pi 
    a2, b2, c2 = 5.0801e-1, 7.7515e-2, 2.2369e-2
    f2 = (a2 * eta ** 2 + b2 * eta + c2) / (m1 + m2) * 2e5 / np.pi 
    a3, b3, c3 = 8.4845e-1, 1.2848e-1, 2.7299e-1
    f3 = (a3 * eta ** 2 + b3 * eta + c3) / (m1 + m2) * 2e5 / np.pi 

    elle = (1. / 2. / np.pi) * f2 / ((ff - f1) ** 2. + f2 ** 2. / 4.)
    w = np.pi * f2 / 2. * (f0 / f1) ** (2. / 3.)

    A = np.sqrt(5.0 / 24.0 / np.pi ** (4.0 / 3.0)) * 3.63e-19 # constants

    mchirpz= mchirp * (1+z)
    #fz = ff / (1+z) # in Hz
    dl = Planck15.luminosity_distance(z).value # in Mpc

    #hc = A * mchirpz ** (5.0 / 6.0) / dl / (1.0 + z) / fz ** (7.0 / 6.0)
    hc = A * mchirpz ** (5.0 / 6.0) / dl / f0 ** (7.0 / 6.0)
    if ff < f0:
        hc *= (ff / f0) ** (-7. / 6.)
    elif ff >= f0 and ff < f1:
        hc *= (ff / f0) ** (-2. / 3.)
    elif ff >= f1 and ff < f3:
        hc *= w * elle
    elif ff > f3:
        hc *= 0.

    return hc



