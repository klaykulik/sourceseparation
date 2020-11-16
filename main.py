import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as cst
from scipy.interpolate import interp1d
from scipy.stats import f as scipy_f

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline

from sourceseparation.telluric_shift import telluric_shift
from sourceseparation.read_hitran import read_hitran, convert


def fit_models(sightlines, iteration, method='leastsq', debug=False, verbose=0):


    if verbose > 1:
        report = True
    elif debug is True:
        report = True
    else:
        report = False

    for i in range(len(sightlines)):
        sightline = sightlines[i]


        print("Fitting sightline...")

        sightline.fit(
            data=sightline.interp_flux,
            x=sightline.grid,
            report=report,
            plot=debug,
            weights=1 / (np.ones_like(sightline.grid) * sightline.sigma),
            method=method
        )

        try:
            sightline.chisqrs.append(sightline.result.chisqr)

        except AttributeError:
            sightline.chisqrs = [sightline.result.chisqr]


        sightline.out = sightline.complete_model.eval(
            data=sightline.interp_flux, params=sightline.result.params, x=sightline.grid
        )

        if debug:
            my_chi_sq = np.sum(((sightline.out - sightline.interp_flux) / sightline.sigma)**2)
            if np.abs(my_chi_sq - sightline.result.chisqr) > 1:
                print()
                print(my_chi_sq)
                print(sightline.result.chisqr)
                raise Exception('Chi-square from lmfit different from calculated Chi-square')
                print()

            # Check other errors
            if 'Could not estimate' in sightline.result.message:
                print()
                raise Exception("Could not estimate error-bars")

            if sightline.result.covar is None:
                print()
                print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
                print(sightline.result.aborted)
                print(sightline.result.message)
                print(sightline.result.errorbars)
                print(sightline.result.covar)
                print(sightline.result.fit_report())
                print(sightline.result.params)
                print()

                raise Exception("No Covariance matrix")


        sightline.telluric_resid = sightline.interp_flux - sightline.out

        if verbose > 0 or debug:
            print('chi-sq: ', sightline.result.chisqr)
            print('reduced chi-sq: ', sightline.result.redchi)

    # %%%%%%%%%%%%%%% Calculate barycentric model and residuals %%%%%%%%%%%%%%%

    for i in range(len(sightlines)):
        sightline = sightlines[i]

        # Shift grid to bary frame
        sightline.bary_shift = (
            sightline.grid
            + (sightline.v_bary / cst.c.to("km/s").value) * sightline.grid
        )

    # find out where all sightlines overlap, and make new grid
    bary_min = np.max([np.min(sightline.bary_shift) for sightline in sightlines])
    bary_max = np.min([np.max(sightline.bary_shift) for sightline in sightlines])

    grid_idx = np.where(
        np.logical_and(sightline.raw_grid > bary_min, sightline.raw_grid < bary_max)
    )
    bary_grid = sightline.raw_grid[grid_idx]

    # interpolate resid and model to spectrum grid
    for i in range(len(sightlines)):
        sightline = sightlines[i]
        sightline.bary_grid = bary_grid

        f = interp1d(sightline.bary_shift, sightline.telluric_resid)
        sightline.bary_resid_interp = f(sightline.bary_grid)

        g = interp1d(sightline.bary_shift, sightline.out)
        sightline.bary_out_interp = g(sightline.bary_grid)


def check_stop(sightlines, iteration, stop_method, freeze_after, verbose=0):
    total_bayes = 0
    total_f = 0

    # #########################
    # BAYES MODEL COMPARISON

    if stop_method == 'bayes':
        posterior_calc(sightlines, verbose)

        for sightline in sightlines:
            if iteration > 0:
                if sightline.posteriors[-1] > 0 and sightline.posteriors[-1] > 0:
                    if sightline.posteriors[-1] / sightline.posteriors[-2] > 1:
                        total_bayes += 1
            else:
                total_bayes += 1

    # #########################
    # F TEST

    elif stop_method == 'f':
        for sightline in sightlines:

            df1 = 4  # number of new parameters in 'better' model
            df2 = len(sightline.grid) - len(sightline.all_pars)

            if iteration > 0:
                num = (sightline.chisqrs[-2] - sightline.chisqrs[-1]) / df1
                denom = sightline.chisqrs[-1] / df2
                f = num / denom
                print("f: ", f)

                alpha = 0.05  # Or whatever you want your alpha to be.
                p_value = scipy_f.cdf(f, df1, df2)
                print("p: ", p_value)
                if p_value > alpha:
                    # Reject the null hypothesis that Var(X) == Var(Y)
                    total_f += 1

    # #########################

    num_to_stop = np.ceil(len(sightlines) / 2.0)

    if total_bayes >= num_to_stop or total_f >= num_to_stop:
        print(
            str(sightlines[0].n_lines)
            + ' lines are better than '
            + str(sightlines[0].n_lines - 1)
        )

        sightline.model_lag_one = sightline.complete_model

        sightline.pars_lag_one = sightline.all_pars





    else:
        print("{} lines does not produce a better model than {} lines".format(
            sightlines[0].n_lines, sightlines[0].n_lines - 1
        ))

        if freeze_after is not None:
            print('Fitting one last time')
            sightline.freeze(unfreeze=True)
            sightline.fit(
                data=sightline.interp_flux,
                params=sightline.pars_lag_one,
                model=sightline.model_lag_one,
                x=sightline.grid,
                report=True,
                plot=True,
                weights=1 / (np.ones_like(sightline.grid) * sightline.sigma),
                method=method
            )

        for sightline in sightlines:
            sightline.separate(data=sightline.interp_flux, x=sightline.grid)
        sys.exit()


def posterior_calc(sightlines, verbose=0):

    # get posterior in each sightline
    for i in range(len(sightlines)):
        sightline = sightlines[i]

        # compute parameter ranges for each line and multiply
        range_terms = []
        for i in range(sightline.n_lines):

            lam_name = sightline.peaks[i].name
            b_name = lam_name.replace('lam_0', 'b')
            d_name = lam_name.replace('lam_0', 'd')
            tau_name = lam_name.replace('lam_0', 'tau_0')
            # print(lam_name, b_name, d_name, tau_name)

            lam_range = (
                sightline.result.params[lam_name].max - sightline.result.params[lam_name].min
            )
            b_range = sightline.result.params[b_name].max - sightline.result.params[b_name].min
            d_range = sightline.result.params[d_name].max - sightline.result.params[d_name].min
            tau_range = (
                sightline.result.params[tau_name].max - sightline.result.params[tau_name].min
            )

            range_term = lam_range * b_range * d_range * tau_range
            range_terms.append(range_term)

        # compute range for continuum parameters
        c_ranges = []
        for name in sightline.all_pars:
            if 'y_' in name:
                c_range = sightline.result.params[name].max - sightline.result.params[name].min
                c_ranges.append(c_range)

        # multiply range term for each line and continuum together
        denominator_term = 1
        for c_range in c_ranges:
            denominator_term = denominator_term * c_range

        for term in range_terms:
            denominator_term = denominator_term * term

        det = 0
        try:
            det = np.linalg.det(sightline.result.covar)
        except np.linalg.LinAlgError:

            if sightline.result.covar is None:
                print()
                print('bad here')
                print(sightline.result.covar)
                print(sightline.result.errorbars)
                print(sightline.result.message)

        a = np.math.factorial(sightline.n_lines) * \
            (4 * np.pi)**((2 * sightline.n_lines + sightline.n_anchors) / 2)

        b = denominator_term * (np.sqrt(det))
        c = np.exp(-sightline.result.chisqr / 2)
        posterior = a / b * c

        try:
            sightline.posteriors.append(posterior)
        except AttributeError:
            sightline.posteriors = [posterior]

        if verbose > 0:
            # print()
            print('a: ', a)
            print('b: ', b)
            print('c: ', c)
            print(sightline.posteriors)


def plot_models(sightlines):
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Fit data
    fig, axs = plt.subplots(5, 2, sharex=True)
    for i in range(len(sightlines)):
        sightline = sightlines[i]


        axs[i, 0].scatter(sightline.grid, sightline.interp_flux,
                          color='k', marker='.', label='Observed data')
        axs[i, 0].plot(sightline.grid, sightline.out, color='r', label='Model')
        axs[i, 0].plot(sightline.grid, sightline.telluric_resid, color='g', label='Residual')
        axs[i, 0].set_ylabel('Flux', fontsize=14)

        axs[i, 1].scatter(sightline.grid, sightline.interp_bary_flux, color='k', marker='.')
        axs[i, 1].plot(sightline.bary_grid, sightline.bary_out_interp, color='r')
        axs[i, 1].plot(sightline.bary_grid, sightline.bary_resid_interp, color='g')

    axs[0, 0].legend()
    axs[4, 0].set_xlabel(r'Wavelength ($\AA$)', fontsize=14)
    axs[4, 1].set_xlabel(r'Wavelength ($\AA$)', fontsize=14)

    axs[0, 0].set_title('Geocentric Reference Frame', fontsize=14)
    axs[0, 1].set_title('Barycentric Reference Frame', fontsize=14)

    # plt.tight_layout()
    plt.show()


def likelihoods(sightlines):
    # %%%%%%%%%%%%%%%%%%%%%%%%% Calculate likelihoods %%%%%%%%%%%%%%%%%%%%%%%%%

    for i in range(len(sightlines)):
        sightline = sightlines[i]

        sightline.log_likelihood = np.log(1 / (sightline.sigma * np.sqrt(2 * np.pi))) + (
            -((sightline.telluric_resid) ** 2 / (2 * sightline.sigma ** 2))
        )

        sightline.bary_log_likelihood = np.log(1 / (sightline.sigma * np.sqrt(2 * np.pi))) + (
            -((sightline.bary_resid_interp) ** 2 / (2 * sightline.sigma ** 2))
        )

    # %%%%%%%%%%%%%%%%%%%%%%%%%% Co-add likelihoods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    coadd = np.ones_like(sightlines[0].grid)
    bary_coadd = np.ones_like(sightlines[0].bary_grid)

    for sightline in sightlines:

        coadd *= sightline.log_likelihood
        bary_coadd *= sightline.bary_log_likelihood


    return coadd, bary_coadd


def plot_likelihoods(sightlines, coadd, bary_coadd):
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Likelihoods
    fig, axs = plt.subplots(1, 2)
    for sightline in sightlines:
        axs[0].plot(sightline.grid, sightline.log_likelihood)
        axs[0].set_xlabel(r'Wavelength ($\AA$)', fontsize=14)


        axs[1].plot(sightline.bary_grid, sightline.bary_log_likelihood)
        axs[1].set_xlabel(r'Wavelength ($\AA$)', fontsize=14)

    # plt.suptitle("Log-Likelihoods in each sightline")
    axs[0].set_title('Geocentric Reference Frame', fontsize=14)
    axs[0].tick_params(axis='both', labelsize=12)
    axs[1].set_title('Barycentric Reference Frame', fontsize=14)
    axs[1].tick_params(axis='both', labelsize=12)

    plt.show()

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # coadded log-likelihoods
    fig, axs = plt.subplots(1, 2)

    axs[0].plot(sightlines[0].grid, coadd)
    axs[0].set_xlabel(r'Wavelength ($\AA$)', fontsize=14)

    axs[1].plot(sightlines[0].bary_grid, bary_coadd)
    axs[1].set_xlabel(r'Wavelength ($\AA$)', fontsize=14)

    axs[0].set_title('Geocentric Reference Frame', fontsize=14)
    axs[0].tick_params(axis='both', labelsize=12)
    axs[1].set_title('Barycentric Reference Frame', fontsize=14)
    axs[1].tick_params(axis='both', labelsize=12)

    plt.show()


def add_next_line(sightlines, coadd, bary_coadd):


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare co-adds %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if np.min(coadd) < np.min(bary_coadd):
        next_source = "Telluric"
    else:
        next_source = "Nontelluric"

    name_next = "line" + str(sightlines[0].n_lines)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add new line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if next_source == "Telluric":
        idx = np.argmin(coadd)

        for sightline in sightlines:
            assert np.min(sightline.interp_flux[idx]) == sightline.interp_flux[idx]

        for i in range(len(sightlines)):
            sightline = sightlines[i]

            # get position and tau_0 of new line
            lam_0 = sightline.grid[idx]
            d = 0.1
            tau_0 = 0.01
            line = {"lam_0": lam_0, "d": d, "tau_0": tau_0}
            sightline.add_line(name=name_next, source=next_source, pars=line)

    elif next_source == "Nontelluric":
        idx = np.argmin(bary_coadd)

        for i in range(len(sightlines)):
            sightline = sightlines[i]

            # get position and tau_0 of new line
            bary_lam_0 = sightline.bary_grid[idx]
            lam_0 = (
                bary_lam_0 - (sightline.v_bary / cst.c.to("km/s").value) * bary_lam_0
            )
            b = 1.2

            # Lifetime broadening of K_I lines - constant value
            K_Gamma = 3.820e7
            K_d = K_Gamma * lam_0**2 / (4 * np.pi * (cst.c.to("cm/s").value * 1e8))


            if sightline.n_lines > 3:
                tau_0 = 0.1
            else:
                tau_0 = 0.05
            line = {"lam_0": lam_0, "b": b, "d": K_d, "tau_0": tau_0}
            sightline.add_line(name=name_next, source=next_source, pars=line)

            sightline.all_pars[next_source + '_' + name_next + '_d'].set(vary=False)


    else:
        print("Something bad happened...")
        sys.exit()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SET PARAMS
method = 'leastsq'
debug = False
verbose = 2
stop_method = 'bayes'
freeze_after = None

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
# TELLURIC WAVELENGTH CORRECTION

zoom_xmin = 7661.5
zoom_xmax = 7669

# Shift spectra
observations = telluric_shift(
    observations, xmin=xmin, xmax=xmax, zoom_xmin=zoom_xmin, zoom_xmax=zoom_xmax, molecule='O2'
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE NEW MODELS

sightlines = []
for i in range(len(observations)):
    sp = observations[i]
    sightline = Sightline(sp, n_anchors=4)
    sightlines.append(sightline)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CALCULATE SIGMA
flat_xmin = 7662.25
flat_xmax = 7663.6

for sightline in sightlines:
    sigma_idx = np.where(
        np.logical_and(sightline.grid > flat_xmin, sightline.grid < flat_xmax)
    )
    sightline.sigma = np.std(sightline.interp_flux[sigma_idx])

    if verbose > 1:
        print('Sigma: ', sightline.sigma)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GET PRIOR DATA

# Re-find O2 lines from HITRAN data
pars_list = convert(read_hitran("telluric_lines_HITRAN.txt", molecule='O2'))

# Set tau_cutoff (0.4, 0.02, or 0.0) (0.4 works, others still under development)
tau_cutoff = 0.4

# Create linelist
linelist = []
for pars in pars_list:
    if (zoom_xmin < pars["lam_0"]) and (pars["lam_0"] < zoom_xmax):
        if pars["tau_0"] > tau_cutoff:
            linelist.append(pars)
linelist.reverse()
if verbose > 1:
    print('Prior Linelist:')
    print(linelist)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ADD PRIOR DATA TO MODELS
for i in range(len(sightlines)):
    sightline = sightlines[i]

    # # Geocentric
    # sightline.add_source(name="Telluric", similar={"b": 0.6})
    for j in range(len(linelist)):
        line = linelist[j].copy()
        name = "line" + str(j)
        line['d'] = 0.1
        sightline.add_line(name=name, source="Telluric", pars=line)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ITERATE

for iteration in range(12):

    print()
    print('################################################################')
    print("Iteration {}, calculating with {} lines".format(iteration, sightlines[0].n_lines))
    print('################################################################')
    print()

    # %%%%%%%%%%%%%%%% Fit geocentric model and get residuals %%%%%%%%%%%%%%%%%
    fit_models(sightlines, iteration=iteration, method=method, debug=debug, verbose=verbose)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check to stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    check_stop(sightlines, iteration=iteration, stop_method=stop_method,
               freeze_after=freeze_after, verbose=1)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_models(sightlines)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Freeze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if freeze_after is not None:
        for sightline in sightlines:
            if sightline.n_lines == freeze_after:
                print("FREEZING sightline after " + str(freeze_after) + " lines")
                sightline.freeze()


    # %%%%%%%%%%%%%%%%%%%% Calculate and coadd likelihoods %%%%%%%%%%%%%%%%%%%%
    coadd, bary_coadd = likelihoods(sightlines)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_likelihoods(sightlines, coadd, bary_coadd)


    # %%%%%%%%%%%%%%%%%%% Compare coaddsd and add new line %%%%%%%%%%%%%%%%%%%%
    add_next_line(sightlines, coadd, bary_coadd)
