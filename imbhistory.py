'''
imbhistory: computing how theoretical merger rates for IMBHs are expected to be observed by different detectors
https://github.com/giacomofragione/imbhistory
'''

import numpy as np
import os, sys

import gwinstr, functions
import detection2 as dt

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
           'qmin' : 1e-2,
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

        red_arr = np.arange(self.z_min, self.z_max, 0.1)
        nev = []
        for i in range(len(red_arr)):
            nev.append(Planck15.comoving_distance(red_arr[i]).value ** 3 / (1. + red_arr[i]) * functions.red_formation(red_arr[i], self.mu, self.sigma))
        nev = np.array(nev)
        ymin, ymax = min(nev), max(nev)

        yfunc, ysample = -1, 0
        while ysample > yfunc:
            ysample = ymin + random.random() * (ymax - ymin)
            xsample = self.z_min + random.random() * (self.z_max - self.z_min)
            yfunc = Planck15.comoving_distance(xsample).value ** 3 / (1. + xsample) * functions.red_formation(xsample, self.mu, self.sigma)

        return xsample


    def compute_snr(self):

        zz = self.ratered()

        m1 = functions.sample_powerlaw(self.mimbh_min, self.mimbh_max, self.slope_mimbh)
        m2 = 0.
        while m2 < 10.:
            qq = functions.sample_powerlaw(self.qmin, self.qmax, self.slope_q)
            m2 = qq * m1

        tlisa = 5. # duration of LISA mission

        ff_min = 0.07 * ((m1 + m2) * (1. + zz) / 100.) ** 0.125 / (m1 * m2 * ((1. + zz) ** 2.) / 100.) ** 0.375 / (tlisa / 4.) ** 0.375

        appr = 'IMRPhenomD'

        psd = 'rob_lisa'
        snr_lisa = (2.0 / np.sqrt(5.)) * np.sqrt(5.) * dt.snr(m1, m2, zz, appr, ff_min, ff_min, psd)
        psd = 'voyager'
        snr_ligov = (2.0 / 5.) * dt.snr(m1, m2, zz, appr, ff_min, ff_min, psd)
        psd = 'et_d'
        snr_et = (2.0 / 5.) * dt.snr(m1, m2, zz, appr, ff_min, ff_min, psd)
        psd = 'ce2'
        snr_ce = (2.0 / 5.) * dt.snr(m1, m2, zz, appr, ff_min, ff_min, psd)

        return(
            zz,
            m1,
            qq,
            snr_lisa,
            snr_ligov,
            snr_et,
            snr_ce,
        )


    def eval(self, nsample):

        zobs_arr, m1obs_arr, qobs_arr, snr_lisa_arr, snr_ligov_arr, snr_et_arr, snr_ce_arr = ([]), ([]), ([]), ([]), ([]), ([]), ([])

        for k in range(nsample):
            zobs, m1obs, qobs, snr_lisa, snr_ligov, snr_et, snr_ce = self.compute_snr()
            zobs_arr = np.append(zobs_arr, zobs)
            m1obs_arr = np.append(m1obs_arr, m1obs)
            qobs_arr = np.append(qobs_arr, qobs)
            snr_lisa_arr = np.append(snr_lisa_arr, snr_lisa)
            snr_ligov_arr = np.append(snr_ligov_arr, snr_ligov)
            snr_et_arr = np.append(snr_et_arr, snr_et)
            snr_ce_arr = np.append(snr_ce_arr, snr_ce)

        return(
            zobs_arr,
            m1obs_arr,
            qobs_arr,
            snr_lisa_arr,
            snr_ligov_arr,
            snr_et_arr,
            snr_ce_arr,
        )


    def __call__(self, nsample):
        ''' Compute mergers detectable by instrument '''

        return self.eval(nsample)


