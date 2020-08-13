import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as cst

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline

from sourceseparation.telluric_shift import telluric_shift
from sourceseparation.read_hitran import read_hitran, convert


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TELLURIC SHIFT

zoom_xmin = 7661.5
zoom_xmax = 7669

# Shift spectra
observations = telluric_shift(
    observations,
    xmin=xmin,
    xmax=xmax,
    zoom_xmin=zoom_xmin,
    zoom_xmax=zoom_xmax
)

for sp in observations:
    # plt.plot(sp.wave, sp.flux)
    plt.plot(sp.grid, sp.interp_flux)
plt.show()

# Ensure all same length
length = [len(sp.interp_flux) for sp in observations]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE NEW MODELS

sightlines = []
bary_sightlines = []
for i in range(len(observations)):
    sp = observations[i]
    sightline = Sightline(sp, n_anchors=4)
    sightlines.append(sightline)

    bary_sightline = Sightline(sp, n_anchors=4)
    bary_sightlines.append(bary_sightline)

# Re-find O2 lines from HITRAN data
pars_list = convert(read_hitran('telluric_lines_HITRAN.txt'))

# Set tau_cutoff (0.4, 0.02, or 0.0)
tau_cutoff = 0.4

# Create linelist
linelist = []
for pars in pars_list:
    if (zoom_xmin < pars['lam_0']) and (pars['lam_0'] < zoom_xmax):
        if pars['tau_0'] > tau_cutoff:
            linelist.append(pars)
linelist.reverse()
print(linelist)

# Add prior telluric lines to models
for i in range(len(sightlines)):
    sightline = sightlines[i]
    bary_sightline = bary_sightlines[i]

    # Geocentric
    sightline.add_source(name='Telluric', similar={'b': 0.6})
    sightline.add_source(name='Barycentric', similar={'b': 0.1})
    for j in range(len(linelist)):
        line = linelist[j].copy()
        name = 'prior_line' + str(j)
        sightline.add_line(name=name, source='Telluric', pars=line)
        par_name_d = 'Telluric_' + name + '_d'
        sightline.model_pars[par_name_d].set(value=0.05)

    # Barycentric
    bary_sightline.add_source(name='Telluric', similar={'b': 0.6})
    bary_sightline.add_source(name='Barycentric', similar={'b': 0.1})
    for j in range(len(linelist)):
        line = linelist[j].copy()
        line['lam_0'] = line['lam_0'] + \
            (bary_sightline.v_bary / cst.c.to("km/s").value) * \
            line['lam_0']
        name = 'prior_line' + str(j)
        bary_sightline.add_line(name=name, source='Telluric', pars=line)
        par_name_d = 'Telluric_' + name + '_d'
        bary_sightline.model_pars[par_name_d].set(value=0.05)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate sigma
flat_xmin = 7667
flat_xmax = 7668

for sightline in sightlines:
    sigma_idx = np.where(np.logical_and(sightline.grid > flat_xmin, sightline.grid < flat_xmax))
    sightline.sigma = np.std(sightline.grid[sigma_idx])
    print('Sigma: ', sightline.sigma)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LOOP

def loop(sightlines, iteration, debug=False):

    # %%%%%%%%%%%%%%%%%%% Fit and get geocentric residuals %%%%%%%%%%%%%%%%%%%%
    if debug:
        fig, axs = plt.subplots(5, 2)

    resids = []
    for i in range(len(sightlines)):
        sightline = sightlines[i]

        print('Fitting geocentric sightline...')
        sightline.fit(data=sightline.interp_flux, x=sightline.grid,
                      report=False, plot=False)
        out = sightline.model.eval(data=sightline.interp_flux,
                                   params=sightline.result.params,
                                   x=sightline.grid)

        resid = out - sightline.interp_flux  # divide or subtract?
        resids.append(resid)

        # Plot
        if debug:
            axs[i, 0].plot(sightline.grid, sightline.interp_flux)
            axs[i, 0].plot(sightline.grid, out)
            axs[i, 0].plot(sightline.grid, resid)

    # %%%%%%%%%%%%%%%%%%% Fit and get barycentric residuals %%%%%%%%%%%%%%%%%%%
    bary_resids = []
    for i in range(len(sightlines)):
        bary_sightline = bary_sightlines[i]

        print('Fitting barycentric sightline...')
        bary_sightline.fit(data=bary_sightline.interp_bary_flux, x=bary_sightline.grid,
                           report=False, plot=False, bary=True)
        bary_out = bary_sightline.model.eval(data=bary_sightline.interp_bary_flux,
                                             params=bary_sightline.bary_result.params,
                                             x=bary_sightline.grid)

        bary_resid = bary_out - bary_sightline.interp_bary_flux  # divide or subtract?
        bary_resids.append(bary_resid)

        # Plot
        if debug:
            axs[i, 1].plot(bary_sightline.grid, bary_sightline.interp_bary_flux)
            axs[i, 1].plot(bary_sightline.grid, bary_out)
            axs[i, 1].plot(bary_sightline.grid, bary_resid)

    if debug:
        plt.show()

    # %%%%%%%%%%%%%%%%%%%%%%%%% Calculate likelihoods %%%%%%%%%%%%%%%%%%%%%%%%%
    likelihoods = []
    bary_likelihoods = []
    for i in range(len(sightlines)):
        resid = resids[i]
        bary_resid = bary_resids[i]
        sightline = sightlines[i]

        log_likelihood = np.log(1 / (sightline.sigma * np.sqrt(2 * np.pi))) + \
            (-((resid) ** 2 / (2 * sightline.sigma ** 2)))

        bary_log_likelihood = np.log(1 / (sightline.sigma * np.sqrt(2 * np.pi))) + \
            (-((bary_resid) ** 2 / (2 * sightline.sigma ** 2)))

        likelihoods.append(log_likelihood)
        bary_likelihoods.append(bary_log_likelihood)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot likelihoods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for likelihood in likelihoods:
        plt.plot(sightline.grid, likelihood)
    plt.title('Geocentric Likelihoods')
    plt.show()
    for bary_likelihood in bary_likelihoods:
        plt.plot(sightlines[0].grid, bary_likelihood)
    plt.title('Barycentric Likelihoods')
    plt.show()

    # %%%%%%%%%%%%%%%%%%%%%%% Co-add residuals and plot %%%%%%%%%%%%%%%%%%%%%%%
    coadd = np.ones_like(sightlines[0].grid)
    bary_coadd = np.ones_like(sightlines[0].grid)

    for i in range(len(sightlines)):
        likelihood = likelihoods[i]
        bary_likelihood = bary_likelihoods[i]

        coadd *= likelihood
        bary_coadd *= bary_likelihood

    fig, axs = plt.subplots(1, 2)

    axs[0].plot(sightlines[0].grid, coadd, label='Geocentric')
    axs[1].plot(sightlines[0].grid, bary_coadd, label='Barycentric')
    plt.suptitle('coadded likelihoods')
    plt.legend()
    plt.show()

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare co-adds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if np.min(coadd) < np.min(bary_coadd):
        next_line = 'Telluric'
    else:
        next_line = 'Barycentric'

    name = 'line' + str(iteration)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add new line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if next_line == 'Telluric':
        idx = np.argmin(coadd)

        for sightline in sightlines:
            assert np.min(sightline.interp_flux[idx]) == sightline.interp_flux[idx]

        for i in range(len(sightlines)):
            sightline = sightlines[i]
            bary_sightline = bary_sightlines[i]

            # get position and tau_0 of new line
            new_lam_0 = sightline.grid[idx]
            tau_0 = 0.2
            line = {'lam_0': new_lam_0, 'd': 0.01, 'tau_0': tau_0}

            # Geocentric
            sightline.add_line(name=name, source=next_line, pars=line)

            # Barycentric
            bary_lam_0 = new_lam_0 + (bary_sightline.v_bary / cst.c.to("km/s").value) * new_lam_0
            bary_line = {'lam_0': bary_lam_0, 'd': 0.01, 'tau_0': tau_0}
            bary_sightline.add_line(name=name, source=next_line, pars=bary_line)
            print('bary_lam_0: ', bary_lam_0)

    elif next_line == 'Barycentric':
        idx = np.argmin(bary_coadd)

        for bary_sightline in bary_sightlines:
            assert np.min(bary_sightline.interp_flux[idx]) == bary_sightline.interp_flux[idx]

        for i in range(len(sightlines)):
            sightline = sightlines[i]
            bary_sightline = bary_sightlines[i]

            # get position and tau_0 of new line
            new_lam_0 = bary_sightline.grid[idx]
            print('new_lam_0: ', new_lam_0)
            tau_0 = 0.2
            bary_line = {'lam_0': new_lam_0, 'd': 0.01, 'tau_0': tau_0}

            # Geocentric
            lam_0 = new_lam_0 - (sightline.v_bary / cst.c.to("km/s").value) * new_lam_0
            line = {'lam_0': lam_0, 'd': 0.01, 'tau_0': tau_0}
            sightline.add_line(name=name, source=next_line, pars=line)
            print('lam_0: ', lam_0)

            # Barycentric
            bary_sightline.add_line(name=name, source=next_line, pars=bary_line)

    else:
        print('Something bad happened...')
        sys.exit()


for i in range(4):
    loop(sightlines, i, debug=True)






































