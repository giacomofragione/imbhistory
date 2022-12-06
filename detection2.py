import pycbc.psd
import pycbc.filter
import pycbc.waveform
import pycbc.types
import astropy.cosmology
import numpy as np
import os
import sys

import random
from datetime import datetime
random.seed(datetime.now())


def snr(m1, m2, z, appr, flow, deltaf, psd, s1x = 0., s1y = 0., s1z = 0., s2x = 0., s2y = 0., s2z = 0.):

    lum_dist = astropy.cosmology.Planck15.luminosity_distance(z).value

    hp, hc = pycbc.waveform.get_fd_waveform(approximant = appr,
                                            mass1 = m1 * (1. + z),
                                            mass2 = m2 * (1. + z),
                                            delta_f = deltaf,
                                            f_lower = flow,
                                            distance = lum_dist,
                    					    spin1x = s1x,
					                        spin1y = s1y,
					                        spin1z = s1z,
					                        spin2x = s2x,
					                        spin2y = s2y,
					                        spin2z = s2z)

    if psd == 'rob_lisa' or psd == 'et_d' or psd == 'voyager' or psd == 'ce2':
            filename = psd+'.txt'
            evaluatedpsd = pycbc.psd.from_txt(filename, len(hp), deltaf, flow, is_asd_file=False)
    else:
            evaluatedpsd = pycbc.psd.analytical.from_string(psd, len(hp), deltaf, flow)

    snr_one = pycbc.filter.matchedfilter.sigma(hp, psd=evaluatedpsd, low_frequency_cutoff=flow)

    return snr_one


def detec_prob_1d(xx):

    if xx > 1.0:
        p = 0.
    else:
        alpha1 = 1.
        a2 = 0.374222
        a4 = 2.04216
        a8 = -2.63948
        p = a2*(1-xx/alpha1)**2.0+a4*(1-xx/alpha1)**4.0+a8 * \
            (1-xx/alpha1)**8.0+(1.0-a2-a4-a8)*(1-xx/alpha1)**10.0

    return p


def detec_prob_3d(xx):

    if xx > 1.4:
        p = 0.
    else:
        alpha1 = 1.4
        a2 = 1.19549
        a4 = 1.61758
        a8 = -4.87024
        p = a2*(1-xx/alpha1)**2.0+a4*(1-xx/alpha1)**4.0+a8 * \
            (1-xx/alpha1)**8.0+(1.0-a2-a4-a8)*(1-xx/alpha1)**10.0

    return p


def star_formation(z):

    return (1+z) ** 2.6 / (1.0 + ((1.0 + z) / 3.2) ** 6.2)


def mean_metallicity(z):

    return 0.153 - 0.074 * z ** 1.34


def star_formation_nsc(z, z_nsc, sigma_nsc):

    return np.exp(-(z - z_nsc) ** 2 / (2.0 * sigma_nsc ** 2))







