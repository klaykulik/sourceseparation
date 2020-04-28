import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from edibles.edibles.utils import EdiblesSpectrum
from edibles.edibles.models import createCont
from edibles.edibles.models import Sightline
from edibles.edibles import fit


def guessNumLines(spectrum):

    cont = createCont((spectrum.wave, spectrum.flux), n_points=3)

    sightline = Sightline(star_name='HD170740', cont=cont)

    fit_m = fit("HD170740", (spectrum.wave, spectrum.flux), sightline.model, silent=True)
    model = fit_m(spectrum.wave)

    resid = np.abs(spectrum.flux - model)

    mean = np.mean(resid)
    std = np.std(spectrum.flux)

    above = resid > (std + mean)

    y_new = np.zeros_like(spectrum.flux)

    y_new[above] = resid[above]
    plt.plot(spectrum.wave, resid)

    plt.hlines(mean + std, xmin=7661.5, xmax=7669.0)

    prominence = (np.max(spectrum.flux) - np.min(spectrum.flux)) * 0.3
    peaks, _ = find_peaks(resid, prominence=prominence)
    peak_wavelengths = [spectrum.wave.iloc[i] for i in peaks]
    peak_fluxes = y_new[peaks]
    # print(peak_wavelengths)
    print(peak_fluxes)

    plt.scatter(peak_wavelengths, peak_fluxes)
    plt.show()







if __name__ == "__main__":

    file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp1 = EdiblesSpectrum(file1)

    xmin = 7661.5
    xmax = 7669.0

    spectrum = sp1.getSpectrum(xmin=xmin, xmax=xmax)

    guessNumLines(spectrum)
