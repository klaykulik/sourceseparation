"""
Version 3 of separation code

"""

import numpy as np
import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum

from wavelength_corr import correctWavelength


FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"
sp1 = EdiblesSpectrum(FILE1)
subset1 = sp1.getSpectrum(xmin=7661.5, xmax=7669)
sigma1 = 0.005
subset1.flux = subset1.flux / np.max(subset1.flux)

sp2 = EdiblesSpectrum(FILE2)
subset2 = sp2.getSpectrum(xmin=7661.5, xmax=7669)
sigma2 = 0.005
subset2.flux = subset2.flux / np.max(subset2.flux)

sp3 = EdiblesSpectrum(FILE3)
subset3 = sp3.getSpectrum(xmin=7661.5, xmax=7669)
sigma3 = 0.005
subset3.flux = subset3.flux / np.max(subset3.flux)

sp4 = EdiblesSpectrum(FILE4)
subset4 = sp4.getSpectrum(xmin=7661.5, xmax=7669)
sigma4 = 0.005
subset4.flux = subset4.flux / np.max(subset4.flux)

sp5 = EdiblesSpectrum(FILE5)
subset5 = sp5.getSpectrum(xmin=7661.5, xmax=7669)
sigma5 = 0.005
subset5.flux = subset5.flux / np.max(subset5.flux)


shifted_wavelength = correctWavelength(
    observations=[sp1, sp2, sp3, sp4, sp5], xmin=7661.5, xmax=7669,
    silent_fit=True, silent_plot=True
)
print(shifted_wavelength)

print(subset1.wave.iloc[45])

shifted1 = list(subset1.wave + shifted_wavelength[0])
shifted2 = list(subset2.wave + shifted_wavelength[0])
shifted3 = list(subset3.wave + shifted_wavelength[0])
shifted4 = list(subset4.wave + shifted_wavelength[0])
shifted5 = list(subset5.wave + shifted_wavelength[0])

print(type(shifted1))
print(shifted1[45])


# for i in range(5):
#     subset1 = (np.add(subset1[0], shifted_wavelength[i]), subset1.flux)
#     subset2 = (np.add(subset2.wave, shifted_wavelength[i]), subset2.flux)
#     subset3 = (np.add(subset3.wave, shifted_wavelength[i]), subset3.flux)
#     subset4 = (np.add(subset4.wave, shifted_wavelength[i]), subset4.flux)
#     subset5 = (np.add(subset5.wave, shifted_wavelength[i]), subset5.flux)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FAKE

print(len(subset1.wave))
print(len(subset2.wave))
print(len(subset3.wave))
print(len(subset4.wave))
print(len(subset5.wave))

print("LOWEST:")
lowest = np.min([
    len(subset1.wave),
    len(subset2.wave),
    len(subset3.wave),
    len(subset4.wave),
    len(subset5.wave)
])

print(lowest)


grid = np.arange(lowest)


fakewave1 = shifted1[0:lowest]
fakewave2 = shifted2[0:lowest]
fakewave3 = shifted3[0:lowest]
fakewave4 = shifted4[0:lowest]
fakewave5 = shifted5[0:lowest]

fakeflux1 = list(subset1.flux[0:lowest])
fakeflux2 = list(subset2.flux[0:lowest])
fakeflux3 = list(subset3.flux[0:lowest])
fakeflux4 = list(subset4.flux[0:lowest])
fakeflux5 = list(subset5.flux[0:lowest])



# Show the spectrum
print(sp1.target)
print("Barycentric Velocity is", sp1.v_bary)
plt.figure()
plt.plot(shifted1, subset1.flux, label="Night1", marker='.')
plt.plot(shifted2, subset2.flux, label="Night2", marker='.')
plt.plot(shifted3, subset3.flux, label="Night3", marker='.')
plt.plot(shifted4, subset4.flux, label="Night4", marker='.')
plt.plot(shifted5, subset5.flux, label="Night5", marker='.')
plt.legend()
plt.title('Real Data')

# plt.show()


realdiff = []
for i in range(lowest):
    realdiff.append(np.std([subset1.flux.iloc[i], subset2.flux.iloc[i], subset3.flux.iloc[i], subset4.flux.iloc[i], subset5.flux.iloc[i]]))

plt.plot(subset1.wave, realdiff)


plt.figure()
plt.plot(grid, fakeflux1, label="Night1", marker='.')
plt.plot(grid, fakeflux2, label="Night2", marker='.')
plt.plot(grid, fakeflux3, label="Night3", marker='.')
plt.plot(grid, fakeflux4, label="Night4", marker='.')
plt.plot(grid, fakeflux5, label="Night5", marker='.')
plt.legend()
plt.title('Fake Data')


fakediff = []
for i in range(lowest):
    fakediff.append(np.std([subset1.flux.iloc[i], subset2.flux.iloc[i], subset3.flux.iloc[i], subset4.flux.iloc[i], subset5.flux.iloc[i]]))

plt.plot(grid, fakediff)

# plt.show()


something1 = [i * j for i, j in zip(fakeflux1, fakediff)]
something2 = [i * j for i, j in zip(fakeflux2, fakediff)]
something3 = [i * j for i, j in zip(fakeflux3, fakediff)]
something4 = [i * j for i, j in zip(fakeflux4, fakediff)]
something5 = [i * j for i, j in zip(fakeflux5, fakediff)]


plt.scatter(grid, something1, label="Night1", marker='.')
plt.scatter(grid, something2, label="Night2", marker='.')
plt.scatter(grid, something3, label="Night3", marker='.')
plt.scatter(grid, something4, label="Night4", marker='.')
plt.scatter(grid, something5, label="Night5", marker='.')
plt.legend()

plt.show()




# # At each wavelength point:
# for i in range(len(sp1.wave)):
#     x = sp1.wave[i]
#     y = sp1.flux[i]

#     # continuum
#     c = 1
#     # interstellar
#     IS = 1
#     # telluric
#     T = 1

#     # model:
#     y = c * IS * T

#     ####################################

#     # PRIOR
#     # continuum: flat between 0.5 and 1.5
#     c_prior_flux = np.linspace(start=0.5, stop=1.5, num=len(sp1.wave))
#     c_prior_prob = np.ones_like(c_prior_flux)
#     c_max = c_prior_flux[np.argmax(c_prior_prob)]

#     # interstellar: flat between 0.0 and 1.0
#     IS_prior_flux = np.linspace(0.0, 1.0, 1000)
#     IS_prior_prob = np.ones_like(IS_prior_flux)
#     IS_max = IS_prior_flux[np.argmax(IS_prior_prob)]

#     # telluric: flat between 0.0 and 1.0
#     T_prior_flux = np.linspace(0.0, 1.0, 1000)
#     T_prior_prob = np.ones_like(T_prior_flux)
#     T_max = T_prior_flux[np.argmax(T_prior_prob)]

#     ####################################

#     # LILKELIHOOD
#     # continuum:
#     c_like_flux = c_prior_flux
#     c_like_prob = 1 / (sigma1 * np.sqrt(2 * np.pi)) \
#         * np.exp(-((x - c_max) ** 2) / 2 * sigma1 ** 2)

#     # interstellar:
#     IS_like_flux = IS_prior_flux
#     IS_like_prob = 1 / (sigma1 * np.sqrt(2 * np.pi)) \
#         * np.exp(-((x - IS_max) ** 2) / 2 * sigma1 ** 2)

#     # telluric
#     T_like_flux = T_prior_flux
#     T_like_prob = 1 / (sigma1 * np.sqrt(2 * np.pi)) \
#         * np.exp(-((x - T_max) ** 2) / 2 * sigma1 ** 2)











# plt.plot(c_prior_flux, c_prior_prob)
# plt.show()
# plt.scatter(c_max, c_like_prob)
# plt.scatter(IS_max, IS_like_prob)
# plt.scatter(T_max, T_like_prob)
# plt.show()
