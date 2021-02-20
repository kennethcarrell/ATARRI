## Pre-whiten data by subtracting the
## primary frequency and next 9 harmonics
## INPUTS:
##    lspg    - a Lomb-Scargle periodogram (from lightkurve)
##    star_lc - a lightcurve (from lightkurve)
## OUTPUTS:
##    prewhitened periodogram
def doPreWhiten(lspg,star_lc):
    # create model
    model_lc = lspg.model(star_lc.time,lspg.frequency_at_max_power)
    for i in range(2,10):
        model_lc = model_lc + lspg.model(star_lc.time,i*lspg.frequency_at_max_power) - 1.0
    # subtract model from lightcurve and 'renormalize'
    diff2 = star_lc - model_lc + 1.0
    # send back the new periodogram
    return diff2.to_periodogram(maximum_frequency=12.0, oversample_factor=50)

import numpy as np

## Find the maximum and minimum phase from the folded lightcurve
## INPUTS:
##    lc - lightcurve (from lightkurve)
##    lspg - Lomb-Scargle periodogram (from lightkurve)
## OUTPUTS:
##    dictionary of max and min values
def get_rise_time(lc, lspg):
    # Get max value for t0 estimate
    t0 = lc.time[np.argmax(lc.flux)]
    # Fold the lightcurve using period from Lomb-Scargle analysis and t0
    folded = lc.fold(period=lspg.period_at_max_power,epoch_time=t0,normalize_phase=True)
    # Get indices of 10 biggest/smallest values
    tenMax = folded.flux.value.argsort()[-10:]
    tenMin = folded.flux.value.argsort()[:10]
    # Find mean values
    maxMean = np.mean(folded.phase.value[tenMax])
    minMean = np.mean(folded.phase.value[tenMin])
    # Return both max and min
    return dict([
        ('max',maxMean),
        ('min',minMean)
    ])


## Estimate uncertainty in the period
## Taken from upsilon package's extract_features.py
## https://github.com/dwkim78/upsilon/blob/master/upsilon/extract_features/extract_features.py
## INPUTS:
##    lspg - Lomb-Scargle periodogram (from lightkurve)
##    fx_width - window width for searching
## OUTPUTS:
##    half-width of peak (at/below average+std)
def get_period_uncertainty(lspg, fx_width=100):
    # get peak index
    jmax = np.argmax(lspg.power)
    sindex = jmax - fx_width
    eindex = jmax + fx_width
    # edge checks
    if(sindex < 0):
        sindex = 0
    if(eindex > len(lspg.power)-1):
        eindex = len(lspg.power)-1
    # grab values in window
    fx_subset = lspg.frequency.value[sindex:eindex]
    fy_subset = lspg.power.value[sindex:eindex]
    # find mean and standard deviation
    fy_mean = np.median(fy_subset)
    fy_std = np.std(fy_subset)
    # find peak again in subset data
    max_index = np.argmax(fy_subset)
    # get values where power is less than average + std
    index = np.where(fy_subset <= fy_mean + fy_std)[0]
    # find left & right edges (full width)
    lindex = index[(index < max_index)]
    if(len(lindex)==0):
        lindex = 0
    else:
        lindex = lindex[-1]
    rindex = index[(index > max_index)]
    if(len(rindex)==0):
        rindex = len(fy_subset)-1
    else:
        rindex = rindex[0]
    # return half-width
    return (1./fx_subset[lindex] - 1./fx_subset[rindex])/2.

#### Fourier Features Functions
## Fourier function
def fourier_series(pars, x, order):
    sum0 = pars[0]
    for i in range(order):
        sum0 += pars[i*2+1]*np.sin(2*np.pi*(i+1)*x) + pars[i*2+2]*np.cos(2*np.pi*(i+1)*x)
    return sum0

## Residual between data and model
def residuals(pars, x, y, order):
    return y - fourier_series(pars, x, order)

from scipy.optimize import leastsq

## Fit data to model and return features
## INPUT:
##    star_lc - lightcurve (from lightkurve)
##    maxper - period of oscillation
## OUTPUT:
##    dictionary of relevant parameters from fit values
def getFeatures(star_lc,maxper):
    # get data in a usable format
    dates   = star_lc.time.value
    flx     = star_lc.flux.value
    flx_err = star_lc.flux_err.value
    # estimate of magnitude
    mag     = -2.5*np.log(flx) + 42
    err     = flx_err/flx*mag

    # fit Fourier Series of order 3
    order = 3
    # initial guess
    p0 = np.ones(order*2 + 1)
    date_period = (dates%maxper)/maxper
    p1, success = leastsq(residuals, p0, args=(date_period, mag, order))
    # from Petersen, J.O., 1986, A&A (via upsilon package)
    return dict(
        amplitude=np.sqrt(p1[1]**2 + p1[2]**2),
        r21=np.sqrt(p1[3]**2 + p1[4]**2)/np.sqrt(p1[1]**2 + p1[2]**2),
        r31=np.sqrt(p1[5]**2 + p1[6]**2)/np.sqrt(p1[1]**2 + p1[2]**2),
        phase=np.arctan(-p1[1] / p1[2]),
        phi21=np.arctan(-p1[3] / p1[4]) - 2.*np.arctan(-p1[1] / p1[2]),
        phi31=np.arctan(-p1[5] / p1[6]) - 3.*np.arctan(-p1[1] / p1[2])
    )
