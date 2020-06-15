import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as cst
import scipy.stats as stats

from edibles.utils.edibles_spectrum import EdiblesSpectrum

from interpolate import interpolate


FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

xmin = 7661.5
xmax = 7669

sp1 = EdiblesSpectrum(FILE1)
subset1 = sp1.getSpectrum(xmin=xmin, xmax=xmax)
sp2 = EdiblesSpectrum(FILE2)
subset2 = sp2.getSpectrum(xmin=xmin, xmax=xmax)
sp3 = EdiblesSpectrum(FILE3)
subset3 = sp3.getSpectrum(xmin=xmin, xmax=xmax)
sp4 = EdiblesSpectrum(FILE4)
subset4 = sp4.getSpectrum(xmin=xmin, xmax=xmax)
sp5 = EdiblesSpectrum(FILE5)
subset5 = sp5.getSpectrum(xmin=xmin, xmax=xmax)

spectra = [sp1, sp2, sp3, sp4, sp5]
iwave, fluxes = interpolate(spectra, xmin=xmin, xmax=xmax)


for i in range(len(fluxes)):
    fluxes[i] = fluxes[i] / np.median(fluxes[i])


# [plt.plot(iwave, flux, marker='.') for flux in fluxes]
# plt.show()


# ######################## first point

# ### MODEL
def model(flux, tell):
    '''
    Args:
        flux (float): The observed flux value
        tell (float): The telluric transmission, 0<tell<1

    Returns:
        cont (float): The flux before passing through the atmosphere
    '''
    cont = tell / flux
    return cont


# ### PRIOR - probability model is correct, or, what we think the initial values should be
fig, axs = plt.subplots(1, 2, sharey=True)
plt.suptitle("Priors")

# Tell should be between 0 and 1 - represents transmission coefficient
# lets say uniform for now?
telluric = np.linspace(0, 1)
tell_prob = stats.uniform.pdf(telluric)
axs[0].plot(telluric, tell_prob)
axs[0].fill_between(telluric, 0, tell_prob, alpha=0.3)
axs[0].set_ylabel("PDF")
axs[0].set_xlabel("Telluric Transmission")

# cont_flux doesnt really have a necessary range - probably near 1
# I DONT THINK YOU NEED THIS...
cont_flux = np.linspace(0, 1)
cont_prob = stats.beta.pdf(cont_flux, 4, 2)
cont_flux += 0.2
axs[1].plot(cont_flux, cont_prob)
axs[1].fill_between(cont_flux, 0, cont_prob, alpha=0.3)
axs[1].set_xlabel("Continuum value")

plt.show()


# ### LIKELIHOOD - probability you get the data, GIVEN THE MODEL IS CORRECT

data = [flux[0] for flux in fluxes]
print(data)
tell = 1

likes = []
[likes.append(model(datum, tell)) for datum in data]
print(likes)

fig, axs = plt.subplots(1, 2, sharey=True)
plt.suptitle("Likelihoods")


telluric = np.linspace(0, 1)
tell_prob = stats.uniform.pdf(telluric)
axs[0].plot(telluric, tell_prob)
axs[0].fill_between(telluric, 0, tell_prob, alpha=0.3)
axs[0].vlines(tell, ymin=0, ymax=1, color='r')
axs[0].set_ylabel("PDF")
axs[0].set_xlabel("Telluric Transmission")


cont_flux = np.linspace(0, 1)
cont_prob = stats.beta.pdf(cont_flux, 4, 2)
cont_flux += 0.2
axs[1].plot(cont_flux, cont_prob)
axs[1].fill_between(cont_flux, 0, cont_prob, alpha=0.3)
axs[1].vlines(likes, ymin=0, ymax=2, color='red')

# axs[1].set_xlabel("Continuum value")

plt.show()

# ### posterior



















# # ##### tell ref frame

# iwaves_bary = []
# for sp in data:
#     v_bary = sp.v_bary

#     iwave_bary = iwave[0] + (sp.v_bary / cst.c.to("km/s").value) * iwave[0]

#     iwaves_bary.append(iwave_bary)


# print(iwaves_bary)


# plt.scatter(iwave[0], 1)
# [plt.scatter(iwave_bary, 2) for iwave_bary in iwaves_bary]
# plt.show()
