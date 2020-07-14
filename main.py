import numpy as np
import matplotlib.pyplot as plt
from specutils.utils.wcs_utils import vac_to_air
from astropy import units as u

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline

from sourceseparation.wavelength_corr import correct

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


for spec in observations:
    plt.plot(spec.wave, spec.flux)
# plt.show()


# observations = correct(observations)


# for spec in observations:
#     plt.plot(spec.wave, spec.flux)
# plt.show()


sightlines = []
for spec in observations:

    sightline = Sightline(spec)
    sightlines.append(sightline)




print(sightlines)


# O_2 data from HITRAN
wavenums = [13041.123638, 13042.947272]
wavelens = [u.AA / (wn / 1e8) for wn in wavenums]
print(wavelens)

air_waves = [vac_to_air(wl, method='Ciddor1996') for wl in wavelens]



print(air_waves)



for air_wave in air_waves:
    plt.scatter(air_wave, 5000)

plt.show()
















