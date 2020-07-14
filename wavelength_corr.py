from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum


def correct(input_spectra):
    '''A function to correct the wavelength shifts


    '''
    # O_2 data from HITRAN
    # wavenums = [13041.123638, 13042.947272, 13050.480755, 13052.322742]
    # wavelens = [1 / (wn / 1e8) for wn in wavenums]


    lefts = []
    for spec in input_spectra:

        prominence = (np.max(spec.flux) - np.min(spec.flux)) * 0.8
        peaks, _ = find_peaks(-spec.flux, prominence=prominence)
        peak_wavelengths = [spec.wave[i] for i in peaks]
        lefts.append(peak_wavelengths[0])

    leftest = np.min(lefts)
    shifts = lefts - leftest


    for i in range(len(input_spectra)):
        spec = input_spectra[i]
        shift = shifts[i]
        spec.wave = spec.wave - shift

    output_spectra = input_spectra

    return output_spectra


if __name__ == "__main__":


    file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"


    xmin = 7661.5
    xmax = 7669.0


    sp1 = EdiblesSpectrum(file1)
    sp1.getSpectrum(xmin, xmax)
    sp2 = EdiblesSpectrum(file2)
    sp2.getSpectrum(xmin, xmax)
    sp3 = EdiblesSpectrum(file3)
    sp3.getSpectrum(xmin, xmax)
    sp4 = EdiblesSpectrum(file4)
    sp4.getSpectrum(xmin, xmax)
    sp5 = EdiblesSpectrum(file5)
    sp5.getSpectrum(xmin, xmax)


    specs = [sp1, sp2, sp3, sp4, sp5]


    for spec in specs:
        plt.plot(spec.wave, spec.flux)
    plt.show()


    new_specs = correct(specs)


    for spec in new_specs:
        plt.plot(spec.wave, spec.flux)
    plt.show()
