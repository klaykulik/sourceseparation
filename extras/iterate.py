import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from specutils.utils.wcs_utils import vac_to_air

from edibles.sightline import Sightline
from edibles.utils.edibles_spectrum import EdiblesSpectrum

from sourceseparation.wavelength_corr import correct
from sourceseparation.interpolate import interpolate

def iterate(sightlines, counter, resid=None):


    # TELLURIC

    t_wave, fluxes = interpolate(sightlines)

    t_coadd = np.ones_like(t_wave)
    for flux in fluxes:
        t_coadd /= flux

    prominence = np.std(t_coadd)
    peaks, _ = find_peaks(t_coadd, prominence=prominence)
    peak_wavelengths = [t_wave[i] for i in peaks]
    peak_fluxes = [t_coadd[i] for i in peaks]

    plt.plot(t_wave, t_coadd)
    plt.hlines(y=np.std(t_coadd), xmin=7662, xmax=7668)
    plt.scatter(peak_wavelengths, peak_fluxes, marker='x')
    plt.show()



    # BARYCENTRIC
    b_wave, fluxes = interpolate(sightlines, bary=True)

    b_coadd = np.ones_like(b_wave)
    for flux in fluxes:
        b_coadd /= flux

    prominence = np.std(b_coadd)
    peaks, _ = find_peaks(b_coadd, prominence=prominence)
    peak_wavelengths = [b_wave[i] for i in peaks]
    peak_fluxes = [b_coadd[i] for i in peaks]

    plt.plot(b_wave, b_coadd)
    plt.hlines(y=np.std(b_coadd), xmin=7662, xmax=7668)
    plt.scatter(peak_wavelengths, peak_fluxes, marker='x')
    plt.show()












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
    observations = [sp1, sp2, sp3, sp4, sp5]



    observations = correct(observations)



    sightlines = []
    for spec in observations:

        sightline = Sightline(spec)

        sightlines.append(sightline)


    print(sightlines[0].v_bary)

    for i in range(1):

        iterate(sightlines, counter=i)
