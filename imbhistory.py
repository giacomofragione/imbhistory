'''
imbhistory: computing how theoretical merger rates for IMBHs are expected to be observed by different detectors
https://github.com/giacomofragione/imbhistory
'''

import numpy as np
import os, sys

import gwinstr, functions

import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

import random
from datetime import datetime
random.seed(datetime.now())

########################################################3

__author__ = "Giacomo Fragione"
__license__ = "MIT"
__version__ = "1.0.1"
__email__ = "giacomo.fragione@northwestern.edu"
this_module='imbhistory'

########################################################3


#Defaults values
defaults={ 'directory' : os.path.dirname(__file__),
           'z_min' : 0,
           'z_max' : 10,
           'mu' : 2,
           'sigma' : 1.5,
           'mimbh_min' : 1e2,
           'mimbh_max' : 1e4,
           'slope_mimbh' : -1,
           'qmin' : 1e-3,
           'qmax' : 1,
           'slope_q' : -1}

class imbhistory(object):
    '''
    Compute cosmological merger rates for supermassive black holes - intermediate-mass black holes binaries

    Usage:
        p=imbhistory()
        zobs, m1obs, qobs, snr = p(nsample)

    Parameters:
        nsample # number of samples
        directory # directory
        z_min # minimum redshift
        z_max # maximum_redshift
        mu # mean redshift distribution
        sigma # dispersion redshift distribution
        mimbh_min # minimum IMBH mass
        mimibh_max # maximum IMBH mass
        slope_mimbh # slope of the power law describing IMBH mass
        qmin # minimum mass ratio
        qmax # maximum mass ratio
        slope_q # slope of the power law describing mass ratio

    Returns:
        zobs # detectable redshift
        m1obs # detectable primary mass
        qobs # detectable mass ratio - m2obs = m1obs * qobs
        snr # signal-to-noise ratio
    '''


    def __init__(self,  directory=defaults['directory'],
                        z_min=defaults['z_min'],
                        z_max=defaults['z_max'],
                        mu=defaults['mu'],
                        sigma=defaults['sigma'],
                        mimbh_min=defaults['mimbh_min'],
                        mimbh_max=defaults['mimbh_max'],
                        slope_mimbh=defaults['slope_mimbh'],
                        qmin=defaults['qmin'],
                        qmax=defaults['qmax'],
                        slope_q=defaults['slope_q']):

        self.directory = directory
        self.z_min = z_min
        self.z_max = z_max
        self.mu = mu
        self.sigma = sigma
        self.mimbh_min = mimbh_min
        self.mimbh_max = mimbh_max
        self.slope_mimbh= slope_mimbh
        self.qmin = qmin
        self.qmax = qmax
        self.slope_q= slope_q


    def ratered(self):

        red = -1
        while red < self.z_min or red > self.z_max: 
            red = np.random.normal(self.mu, self.sigma, 1)[0]

        return red


    def compute_snr(self):

        zz = self.ratered()

        m1 = functions.sample_powerlaw(self.mimbh_min, self.mimbh_max, self.slope_mimbh)
        qq = functions.sample_powerlaw(self.qmin, self.qmax, self.slope_q)
        m2 = qq * m1

        tlisa = 5.

        ff_min = 0.04 * ((m1 + m2) / 100.) ** 0.125 / (m1 * m2 / 100.) ** 0.375 / (tlisa / 4.) ** 0.375
        ff_min = np.log10(ff_min)
        ff_min_lisa = max(ff_min,-5)
        ff_max_lisa = 0
        ff_min_ligo = max(ff_min,1)
        ff_max_ligo = 4
        ff_min_et = max(ff_min,0)
        ff_max_et = 4
        ff_min_decigo = max(ff_min,-3)
        ff_max_decigo = 2
        freq_arr_lisa = np.logspace(ff_min_lisa, ff_max_lisa, 100, endpoint=True)
        freq_arr_ligo = np.logspace(ff_min_ligo, ff_max_ligo, 100, endpoint=True)
        freq_arr_et = np.logspace(ff_min_et, ff_max_et, 100, endpoint=True)
        freq_arr_decigo = np.logspace(ff_min_decigo, ff_max_decigo, 100, endpoint=True)

        sum = 0.0
        freq_arr = frwq_arr_lisa
        for i in range(len(freq_arr)-1):

            freq = (freq_arr[i+1] + freq_arr[i]) / 2.
            deltaf = freq_arr[i+1] - freq_arr[i]
            sum += deltaf * gwinstr.hc(freq, m1, m2, zz) ** 2.0 / gwinstr.lisa_noise(freq)
        snr_lisa = 4.0 / np.sqrt(5.0) * np.sqrt(sum)

        sum = 0.0
        freq_arr = frwq_arr_ligo
        for i in range(len(freq_arr)-1):

            freq = (freq_arr[i+1] + freq_arr[i]) / 2.
            deltaf = freq_arr[i+1] - freq_arr[i]
            sum += deltaf * gwinstr.hc(freq, m1, m2, zz) ** 2.0 / gwinstr.ligo_noise(freq)
        snr_ligo = 4.0 / np.sqrt(5.0) * np.sqrt(sum)

        sum = 0.0
        freq_arr = frwq_arr_et
        for i in range(len(freq_arr)-1):

            freq = (freq_arr[i+1] + freq_arr[i]) / 2.
            deltaf = freq_arr[i+1] - freq_arr[i]
            sum += deltaf * gwinstr.hc(freq, m1, m2, zz) ** 2.0 / gwinstr.et_noise(freq)
        snr_et = 4.0 / np.sqrt(5.0) * np.sqrt(sum)

        sum = 0.0
        freq_arr = frwq_arr_decigo
        for i in range(len(freq_arr)-1):

            freq = (freq_arr[i+1] + freq_arr[i]) / 2.
            deltaf = freq_arr[i+1] - freq_arr[i]
            sum += deltaf * gwinstr.hc(freq, m1, m2, zz) ** 2.0 / gwinstr.decigo_noise(freq)
        snr_decigo = 4.0 / np.sqrt(5.0) * np.sqrt(sum)

        return(
            zz,
            m1,
            qq,
            snr_lisa,
            snr_ligo,
            snr_et,
            snr_decigo,
        )


    def eval(self, nsample):

        zobs_arr, m1obs_arr, qobs_arr, snr_lisa_arr, snr_ligo_arr, snr_et_arr, snr_decigo_arr = ([]), ([]), ([]), ([]), ([]), ([]), ([])

        for k in range(nsample):
            zobs, m1obs, qobs, snr_lisa, snr_ligo, snr_et, snr_decigo = self.compute_snr()
            zobs_arr = np.append(zobs_arr, zobs)
            m1obs_arr = np.append(m1obs_arr, m1obs)
            qobs_arr = np.append(qobs_arr, qobs)
            snr_lisa_arr = np.append(snr_lisa_arr, snr_lisa)
            snr_ligo_arr = np.append(snr_ligo_arr, snr_ligo)
            snr_et_arr = np.append(snr_et_arr, snr_et)
            snr_decigo_arr = np.append(snr_decigo_arr, snr_decigo)

        return(
            zobs_arr,
            m1obs_arr,
            qobs_arr,
            snr_lisa_arr,
            snr_ligo_arr,
            snr_et_arr,
            snr_decigo_arr,
        )


    def __call__(self, nsample):
        ''' Compute mergers detectable by instrument '''

        return self.eval(nsample)


