from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt




def correct(input_spectra):


    for spec in input_spectra:
        # find flat section


        prominence = (np.max(spec.flux) - np.min(spec.flux)) * 0.8
        peaks, _ = find_peaks(-spec.flux, prominence=prominence)
        print(peaks)
        peak_wavelengths = [spec.wave.iloc[i] for i in peaks]
        peak_fluxes = [spec.flux.iloc[i] for i in peaks]



        plt.plot(spec.wave, spec.flux)
        plt.scatter(peak_wavelengths, peak_fluxes)
        plt.hlines(y=np.max(spec.flux) - prominence, xmin=7660, xmax=7670)
        plt.show()








    # return output_spectra





if __name__ == "__main__":



    from edibles.utils.edibles_spectrum import EdiblesSpectrum


    file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"



    xmin = 7661.5
    xmax = 7669.0



    sp1 = EdiblesSpectrum(file1).getSpectrum(xmin, xmax)
    print(type(sp1.wave.iloc[0]))


    correct([sp1])
