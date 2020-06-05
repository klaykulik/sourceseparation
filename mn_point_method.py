import numpy as np
import matplotlib.pyplot as plt
import pymultinest

from edibles.utils.edibles_spectrum import EdiblesSpectrum



FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"


sp1 = EdiblesSpectrum(FILE1)
subset1 = sp1.getSpectrum(xmin=7661.5, xmax=7669)
# subset1.flux = subset1.flux / np.max(subset1.flux)

sp2 = EdiblesSpectrum(FILE2)
subset2 = sp2.getSpectrum(xmin=7661.5, xmax=7669)
# subset2.flux = subset2.flux / np.max(subset2.flux)

sp3 = EdiblesSpectrum(FILE3)
subset3 = sp3.getSpectrum(xmin=7661.5, xmax=7669)
# subset3.flux = subset3.flux / np.max(subset3.flux)

sp4 = EdiblesSpectrum(FILE4)
subset4 = sp4.getSpectrum(xmin=7661.5, xmax=7669)
# subset4.flux = subset4.flux / np.max(subset4.flux)

sp5 = EdiblesSpectrum(FILE5)
subset5 = sp5.getSpectrum(xmin=7661.5, xmax=7669)
# subset5.flux = subset5.flux / np.max(subset5.flux)


def model(flux, T):
    # flux = fstar * T
    return flux / T


# create prior 0<T<1
def prior(cube, ndim, nparams):
    cube[0] = cube[0] * 1.5     # uniform prior between 0.5:1.5
    cube[1] = cube[1]           # uniform prior between 0:1


def loglike(cube, ndim, nparams):
    flux, T = cube[0], cube[1]
    # pos1, width, height1 = cube[0], cube[1], cube[2]

    ymodel = model(flux, T)
    # ymodel = model(pos1, width, height1, height2)

    noise = 0.01
    loglikelihood = (-0.5 * ((ymodel - ydata) / noise)**2).sum()

    return loglikelihood


# analyse the file given as first argument
# datafile = sys.argv[1]
# ydata = numpy.loadtxt(datafile)
ydata = subset1.flux


# analyse with 1 gaussian

# number of dimensions our problem has
parameters = ["flux", "T"]
n_params = len(parameters)

# run MultiNest
pymultinest.run(loglike, prior, n_params, outputfiles_basename='out/', resume = False, verbose = True)