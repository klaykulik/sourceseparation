from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy import constants as cst

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models.create_model import createCont
from edibles.models.model import Sightline
from edibles.fitter import fit

from sourceseparation.wavelength_corr import correctWavelength
from interpolate import interpolate


def coadd(fluxes, xmin, xmax):
    """
    Args:
        observations (list): A list of fluxes
        xmin (float): Minimum wavelength value (Angstroms)
        xmax (float): Maximum wavelength value (Angstroms)

    """
    coadded = np.ones_like(iwave)

    for flux in fluxes:
        coadded *= flux

    return iwave, coadded


def iteration(spectra, fluxes):


    return


if __name__ == "__main__":

    # ############################################################################################
    # INIT
    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

    xmin = 7661.5
    xmax = 7669

    sp1 = EdiblesSpectrum(FILE1)
    sp2 = EdiblesSpectrum(FILE2)
    sp3 = EdiblesSpectrum(FILE3)
    sp4 = EdiblesSpectrum(FILE4)
    sp5 = EdiblesSpectrum(FILE5)
    spectra = [sp1, sp2, sp3, sp4, sp5]

    shifts = correctWavelength(
        spectra, xmin, xmax, silent_fit=True, silent_plot=False
    )

    for i in range(len(spectra)):
        spectra[i].wave = spectra[i].wave + shifts[i]

    # ############################################################################################
    # INTERPOLATE
    iwave, fluxes = interpolate(spectra, xmin=xmin, xmax=xmax)

    # ############################################################################################
    # NORMALIZE
    for i in range(len(fluxes)):
        fluxes[i] = fluxes[i] / np.median(fluxes[i])

    # ############################################################################################
    # PLOT
    [plt.plot(iwave, flux, marker='.') for flux in fluxes]
    plt.show()

    # ############################################################################################
    # COADD SPECTRA

    fluxes = [flux - 1 for flux in fluxes]

    iwave, coadded = coadd(fluxes, xmin, xmax)

    # ############################################################################################
    # FIND PEAKS
    sigma = np.std(coadded)

    prominence = sigma
    peaks, _ = find_peaks(-coadded, prominence=prominence)
    peak_wavelengths = [iwave[i] for i in peaks]
    peak_fluxes = [coadded[i] for i in peaks]

    # ############################################################################################
    # PLOT
    plt.plot(iwave, coadded, label='Coadded spectra', c='k', marker='.')
    plt.hlines(y=-sigma, xmin=xmin, xmax=xmax, label='Sigma', color='b')
    plt.scatter(peak_wavelengths, peak_fluxes, marker='x', label='Peaks', c='r')
    plt.legend()
    plt.show()

    # ############################################################################################
    # FIT

    fluxes = [flux + 1 for flux in fluxes]

    tell_models = []
    for i in range(len(spectra)):

        sp = spectra[i]
        data = iwave, fluxes[i]

        cont = createCont(data, n_points=4)
        sightline = Sightline(star_name=sp.target, cont=cont)
        sightline.addSource(source_name="Telluric", b=0.001, d=0.04)

        for j in range(len(peak_wavelengths)):

            name = "tell_" + str(i)
            sightline.addLine(name=name, lam_0=peak_wavelengths[j], tau_0=0.6)

        fit_model = fit(sp.target, data, sightline.model, breakdown=False, silent=False)

        tell_models.append(fit_model(iwave))

    # ############################################################################################
    # PLOT
    for i in range(len(spectra)):
        plt.subplot(np.floor(len(spectra) / 2), 3, i + 1)

        plt.plot(iwave, fluxes[i])
        plt.plot(iwave, tell_models[i])
        plt.plot(iwave, fluxes[i] - tell_models[i])
        plt.title(spectra[i].date[0:10])

    plt.show()

    # #########################################################################
    # RESIDUALS
    resids = []
    for i in range(len(spectra)):
        resids.append(fluxes[i] - tell_models[i])

    # #########################################################################
    # PLOT
    for i in range(len(spectra)):
        plt.subplot(np.floor(len(spectra) / 2), 3, i + 1)

        plt.plot(iwave, resids[i])
        # plt.plot(iwave, tell_models[i])
        # plt.plot(iwave, fluxes[i] - tell_models[i])
        plt.title(spectra[i].date[0:10])

    plt.show()

    # #########################################################################
    # ITER 2

    # #########################################################################
    # COADD


    iter2_f = []
    [iter2_f.append(resid) for resid in resids]

    iwave, coadded = coadd(iter2_f, xmin, xmax)

    # ############################################################################################
    # FIND PEAKS
    sigma = np.std(coadded)

    prominence = sigma
    peaks, _ = find_peaks(coadded, prominence=prominence)
    peak_wavelengths = [iwave[i] for i in peaks]
    peak_fluxes = [coadded[i] for i in peaks]

    # ############################################################################################
    # PLOT
    [plt.plot(iwave, f, marker='.') for f in iter2_f]

    plt.plot(iwave, coadded, label='Coadded spectra', c='k', marker='.')
    # plt.hlines(y=sigma, xmin=xmin, xmax=xmax, label='Sigma', color='b')
    plt.scatter(peak_wavelengths, peak_fluxes, marker='x', label='Peaks', c='r')
    # plt.legend()
    plt.show()
