'''
imbhistory: computing how theoretical merger rates for IMBHs are expected to be observed by different detectors
https://github.com/giacomofragione/imbhistory
'''

import numpy as np
import os
import sys

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
           'qmin' : 1e-3,
           'qmax' : 1,
           'instr' : lisa}

class imbhistory(object):
    '''
    Compute cosmological merger rates for supermassive black holes - intermediate-mass black holes binaries

    Usage:
        p=imbhistory()
        zobs, m1obs, qobs = p(nsample)

    Parameters:
        nsample # number of samples
        directory # directory
        z_min # minimum redshift
        z_max # maximum_redshift
        mu # mean redshift distribution
        sigma # dispersion redshift distribution
        mimbh_min # minimum IMBH mass
        mimibh_max # maximum IMBH mass
        qmin # minimum mass ratio
        qmax # maximum mass ratio
        instr # detector name

    Returns:
        zobs # detectable redshift
        m1obs # detectable primary mass
        qobs # detectable mass ratio
    '''


    def __init__(self,  directory=defaults['directory'],
                        z_min=defaults['z_min'],
                        z_max=defaults['z_max'],
                        mu=defaults['mu'],
                        sigma=defaults['sigma'],
                        mimbh_min=defaults['mimbh_min'],
                        mimbh_max=defaults['mimbh_max'],
                        qmin=defaults['qmin'],
                        qmax=defaults['qmax'],
                        instr=defaults['instr']):

        self.directory = directory
        self.z_min = z_min
        self.z_max = z_max
        self.mu = mu
        self.sigma = sigma
        self.mimbh_min = mimbh_min
        self.mimbh_max = mimbh_max
        self.qmin = qmin
        self.qmax = qmax
        self.instr = instr


    def ratered(self, z)

        return np.exp(-(z - self.mu) ** 2. / 2. / self.sigma ** 2.)        


    def eval(self, nsample):

        return None


    def __call__(self, nsample):
        ''' Compute mergers detectable by instrument '''

        return self.eval(nsample)


