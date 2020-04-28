import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from edibles.edibles.utils import EdiblesSpectrum
from edibles.edibles.models import createCont
from edibles.edibles.models import Sightline
from edibles.edibles.fitter import fit


def correctWavelength(observations, xmin, xmax, silent_fit=True, silent_plot=True):

    datas = []
    natural = []
    models = []
    # spacing_differences=[] # not currently used
    datas_corrected = []
    shifts = []



    for i in range(len(observations)):

        obs = observations[i]

        # create wavelength data for each observation
        data = obs.getSpectrum(xmin=xmin, xmax=xmax)
        datas.append(data)

        prominence = (np.max(data['flux']) - np.min(data['flux'])) * 0.2
        peaks, _ = find_peaks(-data['flux'], prominence=prominence)
        peak_wavelengths = data.iloc[peaks]

        if silent_fit is False:
            plt.plot(data['wave'], data['flux'])
            plt.scatter(peak_wavelengths['wave'], peak_wavelengths['flux'], marker='x')
            plt.show()

        cont = createCont((data['wave'], data['flux']), n_points=3)
        sightline = Sightline(star_name=obs.target, cont=cont)

        sightline.addSource(source_name="Lines", b=0.01, d=0.06)


        for i in range(len(peak_wavelengths)):
            sightline.addLine(name="line", lam_0=peak_wavelengths.iloc[i,0], tau_0=0.1)
            # sightline.addLine(name="line", lam_0=peak_wavelengths[i], tau_0=0.1)

        model = fit(
            obs.target, (data['wave'], data['flux']), sightline.model, breakdown=False, silent=silent_fit
        )


        models.append(model)
        natural.append(sightline.lines["Lines"][-1].lam_0.val)

        # spacing_differences.append(
        #     np.abs(
        #         sightline.lines["Telluric"][1].lam_0.val
        #         - sightline.lines["Telluric"][0].lam_0.val
        #     )
        # )

    if silent_plot is False:

        fig, axes = plt.subplots(nrows=len(observations), ncols=1, sharex=True)
        for i in range(len(observations)):
            obs = observations[i]
            row = axes[i]
            data = datas[i]

            row.axvline(x=peak_wavelengths.iloc[-1,0], c="r", label=str(peak_wavelengths.iloc[0,0]))
            row.plot(data['wave'], data['flux'])
            ylabel = obs.date[0:10]
            row.set_ylabel(ylabel)

        plt.legend()
        plt.show()

    # print(spacing_differences)

    # shift the wavelength grid
    for i in range(len(observations)):
        obs = observations[i]

        data = datas[i]

        wavelength_shift = np.mean(natural) - natural[i]

        data = (np.add(data['wave'], wavelength_shift), data['flux'])
        datas_corrected.append(data)
        shifts.append(wavelength_shift)

    if silent_plot is False:
        # New plotting
        fig, axes = plt.subplots(nrows=len(observations), ncols=1, sharex=True)
        for i in range(len(observations)):
            obs = observations[i]
            row = axes[i]
            data = datas_corrected[i]

            row.axvline(x=np.mean(natural), c="r", label=str(np.mean(natural)))
            row.axvline(x=7660.42, c="r")

            row.plot(data[0], data[1])
            ylabel = obs.date[0:10]
            row.set_ylabel(ylabel)

        plt.legend()
        plt.show()

    return shifts


if __name__ == "__main__":

    # xmin = 7661.5
    # xmax = 7669
    xmin = 7658
    xmax = 7675

    file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"
    sp1 = EdiblesSpectrum(file1)
    sp2 = EdiblesSpectrum(file2)
    sp3 = EdiblesSpectrum(file3)
    sp4 = EdiblesSpectrum(file4)
    sp5 = EdiblesSpectrum(file5)
    obs = [sp1, sp2, sp3, sp4, sp5]

    datas_corrected = correctWavelength(
        obs, xmin, xmax, silent_fit=True, silent_plot=False
    )

    # file1 = '/HD170740/BLUE_346/HD170740_w346_blue_20140916_O12.fits'
    # file2 = '/HD170740/BLUE_346/HD170740_w346_blue_20150424_O12.fits'
    # file3 = '/HD170740/BLUE_346/HD170740_w346_blue_20160505_O12.fits'
    # file4 = '/HD170740/BLUE_346/HD170740_w346_blue_20160612_O12.fits'
    # file5 = '/HD170740/BLUE_346/HD170740_w346_blue_20170701_O12.fits'
    # sp1 = EdiblesSpectrum(file1)
    # sp2 = EdiblesSpectrum(file2)
    # sp3 = EdiblesSpectrum(file3)
    # sp4 = EdiblesSpectrum(file4)
    # sp5 = EdiblesSpectrum(file5)

    # obs = [sp1, sp2, sp3, sp4, sp5]
    # xmin = 3298
    # xmax = 3305

    # datas_corrected = correctWavelength(obs, xmin, xmax, silent=True)
