import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from edibles.edibles.utils import EdiblesSpectrum
from edibles.edibles.models import createCont
from edibles.edibles.models import Sightline
from edibles.edibles import fit

from sourceseparation.wavelength_corr import correctWavelength


def bayes_separate(observations, xmin, xmax, num_iter):

    for observation in observations:

        obs_subset = observation.getSpectrum(xmin=xmin, xmax=xmax)


        # Telluric

        


        # Barycentric









        plt.plot(obs_subset.wave, obs_subset.flux)
    plt.show()

    return


if __name__ == "__main__":

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # HD170740 KI test region

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

    observations = [sp1, sp2, sp3, sp4, sp5]

    # xmin = 7661.5
    # xmax = 7669.0
    xmin = 7658
    xmax = 7675
    num_iter = 5

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Calculate wavelength shifts

    # shifts = correctWavelength(
    #     raw_observations, xmin, xmax, silent_fit=True, silent_plot=True
    # )

    # for i in range(len(raw_observations)):
    #     raw_observation = raw_observations[i]
    #     shift = shifts[i]

    #     observation = raw_observation

    #     observation.wave = np.add(raw_observation.wave, shift)
    #     observations.append(observation)

    # print("Done adding shifts")

    # for i in range(len(observations)):
    #     observation = observations[i]
    #     raw_observation = raw_observations[i]

    #     plt.plot(raw_observation.wave, raw_observation.flux)
    #     plt.plot(observation.wave, observation.flux)
    #     plt.show()
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sightlines = bayes_separate(observations, xmin=xmin, xmax=xmax, num_iter=num_iter)
