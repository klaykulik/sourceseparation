import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


from edibles.utils.edibles_spectrum import EdiblesSpectrum


def interpolate(data, bary=False):
    ''' Linear interpolation function.

    The purpose of this function is to interpolate multiple spectra onto
    a single wavelength grid. This uses the spacing value CDELT1 from the header of the fits files


    Args:
        data (list): List of EdiblesSpectrum Objects
        bary (bool): If true, will interpolate the barycentric wavelength

    Returns:
        tuple:
            list: Standardized wavelength grid
            list: Each item in this list is  a list of interpolated flux values for that spectrum.
                The length of this list is determined my the number of input spectra.

    '''

    spacing = []
    for spec in data:
        spacing.append(spec.header["CDELT1"])

    _spacing = np.min(spacing)

    if bary:
        xmin = np.max([np.min(spec.bary_wave) for spec in data])
        xmax = np.min([np.max(spec.bary_wave) for spec in data])

        i_wave = np.arange(start=xmin, stop=xmax, step=_spacing)

        fluxes = []
        for spec in data:
            f = interp1d(spec.bary_wave, spec.flux)
            i_flux = f(i_wave)
            fluxes.append(i_flux)


    else:
        xmin = np.max([np.min(spec.wave) for spec in data])
        xmax = np.min([np.max(spec.wave) for spec in data])

        i_wave = np.arange(start=xmin, stop=xmax, step=_spacing)

        fluxes = []
        for spec in data:
            f = interp1d(spec.wave, spec.flux)
            i_flux = f(i_wave)
            fluxes.append(i_flux)

    return i_wave, fluxes


if __name__ == "__main__":

    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

    sp1 = EdiblesSpectrum(FILE1)
    subset1 = sp1.getSpectrum(xmin=7661.5, xmax=7669)

    sp2 = EdiblesSpectrum(FILE2)
    subset2 = sp2.getSpectrum(xmin=7661.5, xmax=7669)

    sp3 = EdiblesSpectrum(FILE3)
    subset3 = sp3.getSpectrum(xmin=7661.5, xmax=7669)

    sp4 = EdiblesSpectrum(FILE4)
    subset4 = sp4.getSpectrum(xmin=7661.5, xmax=7669)

    sp5 = EdiblesSpectrum(FILE5)
    subset5 = sp5.getSpectrum(xmin=7661.5, xmax=7669)

    data = [sp1, sp2, sp3, sp4, sp5]


    i_wave, fluxes = interpolate(data)

    for flux in fluxes:
        plt.plot(i_wave, flux, marker='.')
    plt.show()


    b_wave, fluxes = interpolate(data, bary=True)

    for flux in fluxes:
        plt.plot(b_wave, flux, marker='.')
    plt.show()
