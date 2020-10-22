import matplotlib.pyplot as plt
import numpy as np

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel




from numpy.dual import inv
import numdifftools as ndt




filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"

method = 'least_squares'

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

# for i in range(n_anchors):
#     print(i)
#     if i > 0:

#         pars['x_' + str(i)].set(min=pars['x_' + str(i - 1)].value + 0.001, vary=True)

print(pars)


# pars['x_0'].set(min=0)
# pars['x_1'].set(min=pars['x_0'].value)
# pars['x_2'].set(min=pars['x_1'].value)
# pars['x_3'].set(min=pars['x_2'].value)


# for i in range(n_anchors):
    # if i > 0:
    #     pars['x_%i' % (i)].set(min=pars['x_%i' % (i - 1)].value)



result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
out = model.eval(data=sp.flux, params=result.params, x=sp.wave)
resid = sp.flux - out

# print(result.fit_report())
result.plot_fit()
plt.show()



print(result.covar)
print()
print()
print()
print()


# #################################################################################


voigt1 = VoigtModel(prefix='v1_')
voigt1_pars = voigt1.guess(sp.flux, x=sp.wave)


voigt2 = VoigtModel(prefix='v2_')
voigt2_pars = voigt2.guess(sp.flux, x=sp.wave)


voigt3 = VoigtModel(prefix='v3_')
voigt3_pars = voigt3.guess(sp.flux, x=sp.wave)


model = voigt1 * voigt2 * voigt3
pars = voigt1_pars + voigt2_pars + voigt3_pars


pars['v1_lam_0'].set(value=7664.794, min=3000, max=7665.25)
pars['v1_b'].set(value=0.8, min=0.5, max=5)
pars['v1_d'].set(value=0.04, min=0, max=10)
pars['v1_tau_0'].set(value=0.76, min=0, max=10)


pars['v2_lam_0'].set(value=7665.25, min=7665, max=7665.5)
pars['v2_b'].set(value=1.9, min=0.49, max=5)
pars['v2_d'].set(value=0.002, min=0, max=10)
pars['v2_tau_0'].set(value=0.15, min=0, max=10)


pars['v3_lam_0'].set(value=7665.9, min=7665.5, max=7667)
pars['v3_b'].set(value=0.79, min=0.51, max=5)
pars['v3_d'].set(value=0.04, min=0, max=10)
pars['v3_tau_0'].set(value=0.72, min=0, max=10)





# result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
# out = model.eval(data=sp.flux, params=result.params, x=sp.wave)
# resid = sp.flux - out


# result.plot_fit()
# plt.show()

# print(result.covar)
# print()
# print()
# print()
# print()

# #################################################################################

# print(np.all(np.isfinite(out)))
# print(model.param_names)
# print(result.x)
# print()
# print(result.penalty)
# print()



model = cont_model * voigt1 * voigt2 * voigt3

pars = cont_pars + pars



result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
result.plot_fit()
plt.show()
print("the real answer:")
print(result.covar)

print()
print()
print()
print()


# print(result.aborted)
# print(len(result.residual) > len(result.var_names))




# # print(result.penalty)



# # print(result.params)
# print(result.x)

# # test = np.array([
# #     7.66479465e+03, 8.34800324e-01, 4.35296180e-02, 7.28037043e-01,
# #     7.66527454e+03, 1.95507309e+00, 2.31807517e-03, 1.57555907e-01,
# #     7.66586677e+03, 7.93483416e-01, 4.36312498e-02, 7.63314901e-01
# # ])
# # print(test)


# print()
# print()
# print()





# print(np.all(np.isfinite(out)))
# print(model.param_names)
# print(result.x)
# print()
# print(result.penalty)
# print()





# Hfun = ndt.Hessian(result.penalty, step=1e-4)



# hessian_ndt = Hfun(result.x)

# print(hessian_ndt)

# cov_x = inv(hessian_ndt) * 2.0


# print(len(result.residual))
# print(len(result.var_names))
# print(result.nfev)


# print(type(result.covar))
# print(result.covar)

# # plt.plot(sp.wave, sp.flux)
# # plt.plot(sp.wave, out)
# # plt.plot(sp.wave, resid)
# # plt.show()

#  # if self.nan_policy == 'raise' and not np.all(np.isfinite(model))