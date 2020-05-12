"""
Version 3 of separation code

"""

import numpy as np
import matplotlib.pyplot as plt

from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum


FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"
sp1 = EdiblesSpectrum(FILE1)
sigma1 = 0.005
sp2 = EdiblesSpectrum(FILE2)
sigma2 = 0.005
sp3 = EdiblesSpectrum(FILE3)
sigma3 = 0.005
sp4 = EdiblesSpectrum(FILE4)
sigma4 = 0.005
sp5 = EdiblesSpectrum(FILE5)
sigma5 = 0.005

print(len(sp1.wave))
print(len(sp2.wave))
print(len(sp3.wave))
print(len(sp4.wave))
print(len(sp5.wave))

# Show the spectrum
print(sp1.target)
print("Barycentric Velocity is", sp1.v_bary)
plt.figure()
plt.plot(sp1.wave, sp1.flux, label="Geocentric")
plt.show()

# At each wavelength point:
for i in range(len(sp1.wave)):
    x = sp1.wave[i]
    y = sp1.flux[i]

    # continuum
    c = 1
    # interstellar
    IS = 1
    # telluric
    T = 1

    # model:
    y = c * IS * T

    ####################################

    # PRIOR
    # continuum: flat between 0.5 and 1.5
    c_prior_flux = np.linspace(start=0.5, stop=1.5, num=len(sp1.wave))
    c_prior_prob = np.ones_like(c_prior_flux)
    c_max = c_prior_flux[np.argmax(c_prior_prob)]

    # interstellar: flat between 0.0 and 1.0
    IS_prior_flux = np.linspace(0.0, 1.0, 1000)
    IS_prior_prob = np.ones_like(IS_prior_flux)
    IS_max = IS_prior_flux[np.argmax(IS_prior_prob)]

    # telluric: flat between 0.0 and 1.0
    T_prior_flux = np.linspace(0.0, 1.0, 1000)
    T_prior_prob = np.ones_like(T_prior_flux)
    T_max = T_prior_flux[np.argmax(T_prior_prob)]

    ####################################

    # LILKELIHOOD
    # continuum:
    c_like_flux = c_prior_flux
    c_like_prob = 1 / (sigma1 * np.sqrt(2 * np.pi)) \
        * np.exp(-((x - c_max) ** 2) / 2 * sigma1 ** 2)

    # interstellar:
    IS_like_flux = IS_prior_flux
    IS_like_prob = 1 / (sigma1 * np.sqrt(2 * np.pi)) \
        * np.exp(-((x - IS_max) ** 2) / 2 * sigma1 ** 2)

    # telluric
    T_like_flux = T_prior_flux
    T_like_prob = 1 / (sigma1 * np.sqrt(2 * np.pi)) \
        * np.exp(-((x - T_max) ** 2) / 2 * sigma1 ** 2)











# plt.plot(c_prior_flux, c_prior_prob)
# plt.show()
# plt.scatter(c_max, c_like_prob)
# plt.scatter(IS_max, IS_like_prob)
# plt.scatter(T_max, T_like_prob)
# plt.show()
