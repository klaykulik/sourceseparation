import numpy as np
import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum


file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

xmin = 7661.5
xmax = 7669.0



sp1 = EdiblesSpectrum(file1).getSpectrum(xmin, xmax)
sp2 = EdiblesSpectrum(file2).getSpectrum(xmin, xmax)
sp3 = EdiblesSpectrum(file3).getSpectrum(xmin, xmax)
sp4 = EdiblesSpectrum(file4).getSpectrum(xmin, xmax)
sp5 = EdiblesSpectrum(file5).getSpectrum(xmin, xmax)


observations = [sp1, sp2, sp3, sp4, sp5]


plt.plot(sp1.wave, sp1.flux)
plt.plot(sp2.wave, sp2.flux)
plt.plot(sp3.wave, sp3.flux)
plt.plot(sp4.wave, sp4.flux)
plt.plot(sp5.wave, sp5.flux)

plt.show()






















