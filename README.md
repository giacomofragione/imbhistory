# imbhistory

'''
imbhistory: computing how theoretical merger rates for IMBHs are expected to be observed by different detectors
https://github.com/giacomofragione/imbhistory
'''

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

