import matplotlib.pyplot as plt
import numpy as np
from numpy.dual import inv
import numdifftools as ndt

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel
from edibles.sightline import Sightline

from sourceseparation.read_hitran import read_hitran, convert






filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"


sp = EdiblesSpectrum(filename)
print(sp.target)
sp.getSpectrum(xmin=7661, xmax=7670)
# sp.flux = sp.flux / np.max(sp.flux)

# #################################################################################

n_anchors = 4

cont_model = ContinuumModel(n_anchors=n_anchors)
cont_pars = cont_model.guess(sp.flux, x=sp.wave)
model = cont_model
pars = cont_pars






# result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)




# print(result.covar)
# print()
# print()
# print()



voigt1 = VoigtModel(prefix='v1_')
voigt1_pars = voigt1.guess(sp.flux, x=sp.wave)

voigt2 = VoigtModel(prefix='v2_')
voigt2_pars = voigt2.guess(sp.flux, x=sp.wave)

voigt3 = VoigtModel(prefix='v3_')
voigt3_pars = voigt3.guess(sp.flux, x=sp.wave)

voigt4 = VoigtModel(prefix='v4_')
voigt4_pars = voigt4.guess(sp.flux, x=sp.wave)

voigt5 = VoigtModel(prefix='v5_')
voigt5_pars = voigt5.guess(sp.flux, x=sp.wave)

voigt6 = VoigtModel(prefix='v6_')
voigt6_pars = voigt6.guess(sp.flux, x=sp.wave)

voigt7 = VoigtModel(prefix='v7_')
voigt7_pars = voigt7.guess(sp.flux, x=sp.wave)

voigt8 = VoigtModel(prefix='v8_')
voigt8_pars = voigt8.guess(sp.flux, x=sp.wave)


model = cont_model * voigt1 * voigt2 * voigt3
pars = cont_pars + voigt1_pars + voigt2_pars + voigt3_pars


pars['v1_lam_0'].set(value=7664.794, max=7665.25)
pars['v1_b'].set(value=0.8, min=0.5, max=5)
pars['v1_d'].set(value=0.04, min=0, max=10)
pars['v1_tau_0'].set(value=0.76, min=0, max=10)

pars['v2_lam_0'].set(value=7665.25, min=7665, max=7665.5)
pars['v2_b'].set(value=1.9, min=0.49, max=5)
pars['v2_d'].set(value=0.002, min=0, max=10)
pars['v2_tau_0'].set(value=0.15, min=0, max=10)

pars['v3_lam_0'].set(value=7665.9, min=7665.5)
pars['v3_b'].set(value=0.79, min=0.51, max=5)
pars['v3_d'].set(value=0.04, min=0, max=10)
pars['v3_tau_0'].set(value=0.72, min=0, max=10)






result = model.fit(data=sp.flux, params=pars, x=sp.wave)
result.plot_fit()
plt.show()
print("the real answer:")
print(result.covar)

print()
print()
print()
print()



# # #################################################################################
# # #################################################################################
# # #################################################################################




# zoom_xmin = 7661.5
# zoom_xmax = 7669



# # Re-find O2 lines from HITRAN data
# pars_list = convert(read_hitran('sourceseparation/telluric_lines_HITRAN.txt'))

# # Set tau_cutoff (0.4, 0.02, or 0.0)
# tau_cutoff = 0.4

# # Create linelist
# linelist = []
# for pars in pars_list:
#     if (zoom_xmin < pars['lam_0']) and (pars['lam_0'] < zoom_xmax):
#         if pars['tau_0'] > tau_cutoff:
#             linelist.append(pars)
# linelist.reverse()
# print(linelist)








# sightline = Sightline(sp, n_anchors=4)
# sightline.add_source(name='Telluric', similar={'b': 0.6})
# sightline.add_source(name='Barycentric', similar={'b': 0.1})



# sightline.add_line(name='line1', source='Telluric', pars=linelist[0])
# sightline.model_pars["Telluric_line1_d"].set(value=0.05)
# # sightline.add_line(name='line2', source='Telluric', pars=linelist[1])
# # sightline.model_pars["Telluric_line2_d"].set(value=0.05)


# sightline.fit(data=sightline.interp_flux, x=sightline.grid, report=True, plot=True)


# print(sightline.result.covar)




# FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
# xmin = 7661.5
# xmax = 7669

# sp1 = EdiblesSpectrum(FILE1)
# sp1.getSpectrum(xmin=7661, xmax=7670)

# sightline = Sightline(sp1)







# sightline.fit(report=True, plot=True, method='leastsq')
# print(sightline.result.covar)


















# # Add source
# sightline.add_source('telluric')

# # Add line with auto-guessed params
# sightline.add_line(name='line1', source='telluric')





# sightline.fit(report=True, plot=True, method='leastsq')
# out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
# resid = sp1.flux - out
# print(sightline.result.covar)






# # Add line with auto-guessed params
# sightline.add_line(name='line2', source='telluric', guess_data=resid)

# sightline.model_pars['telluric_line2_lam_0'].set(max=sightline.model_pars['telluric_line1_lam_0'].value - 0.2)
# print(sightline.model_pars)


# sightline.fit(report=True, plot=True, method='leastsq')
# print(sightline.result.covar)
















# # Add line with user defined params
# d = {'d': 0.01, 'tau_0': 0.6, 'lam_0': 7664.8}
# sightline.add_line(name='line2', pars=d, source='telluric')

# # Add line with different source
# # d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7665.2}
# sightline.add_source('interstellar', similar={'b': 1.9})
# sightline.add_line(name='line3', source='interstellar')

# # Add line with no source & user defined pars
# d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7662}
# sightline.add_line(name='line4', pars=d)

# # ###############################################################
# # Fit and plot
# sightline.fit(report=True, plot=True, method='leastsq')

# out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
# resid = sp1.flux - out

# plt.plot(sp1.wave, sp1.flux)
# plt.plot(sp1.wave, out)
# plt.plot(sp1.wave, resid)
# plt.show()
# # ###############################################################

# # Add line using guess_pars, and link parameters together
# sightline.add_line(name='line5', source='interstellar', guess_data=resid)
# sightline.model_pars['interstellar_line5_lam_0'].set(expr='interstellar_line3_lam_0 + 0.091')

# # ###############################################################
# # Fit and plot
# sightline.fit(report=True, plot=True, method='leastsq')

# out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
# resid = sp1.flux - out

# plt.plot(sp1.wave, sp1.flux)
# plt.plot(sp1.wave, out)
# plt.plot(sp1.wave, resid)
# plt.show()
# # ###############################################################



