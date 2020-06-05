import pymultinest

def prior(cube, ndim, nparams):
    cube[0] = cube[0] * 2

def loglikelihood(cube, ndim, nparams):
    return -0.5 * ((cube[0] - 0.2) / 0.1)**2

pymultinest.run(loglikelihood, prior, n_params=1, n_dims=1, outputfiles_basename='out/')