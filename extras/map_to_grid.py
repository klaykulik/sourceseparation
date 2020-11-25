from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import numpy as np

from edibles.utils.functions import make_grid
from edibles.utils.edibles_spectrum import EdiblesSpectrum


def mapToGrid(spec, xmin, xmax, resolution=50000):
    '''
    This function will interpolate a given input spectrum onto
    predetermined grid points.

    '''

    grid = make_grid(lambda_start=xmin, lambda_end=xmax, resolution=resolution)

    spline = CubicSpline(spec.wave, spec.flux)

    new_spec = (grid, spline(grid))

    return new_spec


if __name__ == "__main__":

    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

    xmin = 7661.5
    xmax = 7669

    sp1 = EdiblesSpectrum(FILE1)
    subset1 = sp1.getSpectrum(xmin=xmin, xmax=xmax)
    subset1.flux = subset1.flux / np.max(subset1.flux)

    sp2 = EdiblesSpectrum(FILE2)
    subset2 = sp2.getSpectrum(xmin=xmin, xmax=xmax)
    subset2.flux = subset2.flux / np.max(subset2.flux)

    sp3 = EdiblesSpectrum(FILE3)
    subset3 = sp3.getSpectrum(xmin=xmin, xmax=xmax)
    subset3.flux = subset3.flux / np.max(subset3.flux)

    sp4 = EdiblesSpectrum(FILE4)
    subset4 = sp4.getSpectrum(xmin=xmin, xmax=xmax)
    subset4.flux = subset4.flux / np.max(subset4.flux)

    sp5 = EdiblesSpectrum(FILE5)
    subset5 = sp5.getSpectrum(xmin=xmin, xmax=xmax)
    subset5.flux = subset5.flux / np.max(subset5.flux)

    new_spec1 = mapToGrid(spec=subset1, xmin=xmin, xmax=xmax)
    new_spec2 = mapToGrid(spec=subset2, xmin=xmin, xmax=xmax)
    new_spec3 = mapToGrid(spec=subset3, xmin=xmin, xmax=xmax)
    new_spec4 = mapToGrid(spec=subset4, xmin=xmin, xmax=xmax)
    new_spec5 = mapToGrid(spec=subset5, xmin=xmin, xmax=xmax)

    plt.scatter(subset1.wave, subset1.flux)
    plt.scatter(subset2.wave, subset2.flux)
    plt.scatter(subset3.wave, subset3.flux)
    plt.scatter(subset4.wave, subset4.flux)
    plt.scatter(subset5.wave, subset5.flux)

    plt.plot(new_spec1[0], new_spec1[1], marker='.')
    plt.plot(new_spec2[0], new_spec2[1], marker='.')
    plt.plot(new_spec3[0], new_spec3[1], marker='.')
    plt.plot(new_spec4[0], new_spec4[1], marker='.')
    plt.plot(new_spec5[0], new_spec5[1], marker='.')

    plt.show()
