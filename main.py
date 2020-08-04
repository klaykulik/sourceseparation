import os
import numpy as np
import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum

from sourceseparation.telluric_shift import telluric_shift


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# IMPORT DATA

file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

xmin = 7640
xmax = 7680

sp1 = EdiblesSpectrum(file1)
sp2 = EdiblesSpectrum(file2)
sp3 = EdiblesSpectrum(file3)
sp4 = EdiblesSpectrum(file4)
sp5 = EdiblesSpectrum(file5)
observations = [sp1, sp2, sp3, sp4, sp5]


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CORRECT TELLURIC SHIFT


# for sp in observations:
    # print(len(sp.wave))
    # print(len(sp.bary_wave))
    # print()




observations = telluric_shift(observations, xmin=xmin, xmax=xmax, zoom_min=7661, zoom_max=7670)

for sp in observations:
    # plt.plot(sp.wave, sp.flux)
    plt.plot(sp.grid, sp.interp_flux)
plt.show()





coadd = np.ones_like(observations[0].interp_flux)
for sp in observations:

    coadd /= sp.interp_flux


plt.plot(sp.grid, coadd)
plt.show()
