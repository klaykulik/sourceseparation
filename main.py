import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as cst
from scipy.interpolate import interp1d


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
    observations, xmin=xmin, xmax=xmax, zoom_xmin=zoom_xmin, zoom_xmax=zoom_xmax
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NORMALIZE

# for obs in observations:
#     obs.raw_flux = obs.raw_flux / np.max(obs.raw_flux)
#     obs.flux = obs.flux / np.max(obs.flux)
#     obs.interp_flux = obs.interp_flux / np.max(obs.interp_flux)
#     obs.interp_bary_flux = obs.interp_bary_flux / np.max(obs.interp_bary_flux)


# for obs in observations:
#     plt.plot(obs.wave, obs.flux)
#     plt.plot(obs.grid, obs.interp_flux)
# plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE NEW MODELS

sightlines = []
for i in range(len(observations)):
    sp = observations[i]
    sightline = Sightline(sp, n_anchors=5)
    sightlines.append(sightline)

# Re-find O2 lines from HITRAN data
pars_list = convert(read_hitran("telluric_lines_HITRAN.txt"))

# Set tau_cutoff (0.4, 0.02, or 0.0)
tau_cutoff = 0.4

# Create linelist
linelist = []
for pars in pars_list:
    if (zoom_xmin < pars["lam_0"]) and (pars["lam_0"] < zoom_xmax):
        if pars["tau_0"] > tau_cutoff:
            linelist.append(pars)
linelist.reverse()
print(linelist)

# Add prior telluric lines to models
for i in range(len(sightlines)):
    sightline = sightlines[i]

    # Geocentric
    sightline.add_source(name="Telluric", similar={"b": 0.6})
    # sightline.add_source(name='Barycentric', similar={'b': 1.5})
    for j in range(len(linelist)):
        line = linelist[j].copy()
        name = "prior_line" + str(j)
        sightline.add_line(name=name, source="Telluric", pars=line)
        par_name_d = "Telluric_" + name + "_d"
        sightline.model_pars[par_name_d].set(value=0.05, min=0, max=10)

    # sightline.fit(data=sightline.interp_flux, x=sightline.grid, report=True, plot=True)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate sigma
flat_xmin = 7662.25
flat_xmax = 7663.6

for sightline in sightlines:
    sigma_idx = np.where(
        np.logical_and(sightline.grid > flat_xmin, sightline.grid < flat_xmax)
    )
    sightline.sigma = np.std(sightline.interp_flux[sigma_idx])

    print('Sigma: ', sightline.sigma)
    # print('Weight: ', 1 / (sightline.sigma**2))


    # plt.plot(sightline.grid, sightline.interp_flux)
    # plt.hlines(np.mean(sightline.interp_flux[sigma_idx]), flat_xmin, flat_xmax)
    # plt.hlines(np.mean(sightline.interp_flux[sigma_idx]) + sightline.sigma, flat_xmin, flat_xmax)
    # plt.hlines(np.mean(sightline.interp_flux[sigma_idx]) - sightline.sigma, flat_xmin, flat_xmax)
    # plt.show()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LOOP

def loop(sightlines, iteration, method='leastsq', debug=False, plot=False):

    # %%%%%%%%%%%%%%%% Fit geocentric model and get residuals %%%%%%%%%%%%%%%%%

    telluric_resids = []
    outs = []

    print()
    print()
    print('################################################################')
    print("Iteration ", iteration)
    print('################################################################')



    for i in range(len(sightlines)):
        sightline = sightlines[i]

        print("Fitting sightline...")

        try:

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            fudge_factor = 1
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

            sightline.fit(
                data=sightline.interp_flux,
                x=sightline.grid,
                report=debug,
                plot=debug,
                weights=1 / (np.ones_like(sightline.grid) * sightline.sigma * fudge_factor),
                method=method
            )

            # Make sure chi-square calculated from lmfit is correct
            out = sightline.model.eval(
                data=sightline.interp_flux, params=sightline.result.params, x=sightline.grid
            )
            my_chi_sq = np.sum(((out - sightline.interp_flux) / sightline.sigma / fudge_factor)**2)
            if np.abs(my_chi_sq - sightline.result.chisqr) > 1:
                print()
                raise Exception('Chi-square from lmfit different from calculated Chi-square')
                print()

            # Check other errors
            if 'Could not estimate' in sightline.result.message:
                print()
                raise Exception("Could not estimate error-bars")
                print()

            if sightline.result.covar is None:
                print()
                raise Exception("No Covariance matrix")
                print()

        except Exception:
            print()
            print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            print(sightline.result.aborted)
            print(sightline.result.message)
            print(sightline.result.errorbars)
            print(sightline.result.covar)
            print(sightline.result.params)


            if method == 'leastsq':
                print()
                print("leastsq path")

                from scipy.optimize import leastsq as scipy_leastsq

                lskws = sightline.result.call_kws

                print(lskws)
                print(sightline.result.residual)
                print(sightline.result.init_vals)

                lsout = scipy_leastsq(
                    sightline.result._Minimizer__residual,
                    sightline.result.init_vals,
                    **lskws
                )
                _best, _cov, _infodict, errmsg, ier = lsout

                print(_cov)
                print(_best)

                # print(lsout)
                print()


            elif method == 'least_squares':
                print()
                print("least_squares path")

                from numpy.dual import inv
                import numdifftools as ndt

                Hfun = ndt.Hessian(sightline.result.penalty, step=1e-4)
                hessian_ndt = Hfun(sightline.result.x)

                print(hessian_ndt)

                cov_x = inv(hessian_ndt) * 2.0

                print(cov_x)
                print()


            import traceback
            traceback.print_exc()

            sys.exit()

        out = sightline.model.eval(
            data=sightline.interp_flux, params=sightline.result.params, x=sightline.grid
        )
        outs.append(out)

        telluric_resid = sightline.interp_flux - out
        telluric_resids.append(telluric_resid)

        print('chi-sq: ', sightline.result.chisqr)

    # %%%%%%%%%%%%%%% Calculate barycentric model and residuals %%%%%%%%%%%%%%%
    bary_resids = []
    bary_outs = []
    bary_shifts = []

    for i in range(len(sightlines)):
        sightline = sightlines[i]
        telluric_resid = telluric_resids[i]
        out = outs[i]

        # Shift grid to bary frame
        bary_shift = (
            sightline.grid
            + (sightline.v_bary / cst.c.to("km/s").value) * sightline.grid
        )
        bary_shifts.append(bary_shift)

    # find out where all sightlines overlap, and make new grid
    bary_min = np.max([np.min(bary_shift) for bary_shift in bary_shifts])
    bary_max = np.min([np.max(bary_shift) for bary_shift in bary_shifts])

    grid_idx = np.where(
        np.logical_and(sightline.raw_grid > bary_min, sightline.raw_grid < bary_max)
    )
    bary_grid = sightline.raw_grid[grid_idx]

    # interpolate resid and model to spectrum grid
    for i in range(len(sightlines)):
        sightline = sightlines[i]
        telluric_resid = telluric_resids[i]
        out = outs[i]
        bary_shift = bary_shifts[i]

        f = interp1d(bary_shift, telluric_resid)
        bary_resid_interp = f(bary_grid)
        bary_resids.append(bary_resid_interp)

        g = interp1d(bary_shift, out)
        bary_out_interp = g(bary_grid)
        bary_outs.append(bary_out_interp)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check to stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stop(sightlines)

    total = 0
    for sightline in sightlines:

        if len(sightline.peaks) > 4:
            if sightline.answers[-1] / sightline.answers[-2] > 1:
                total += 1
        else:
            total += 1

    if total >= 3:
        print(
            str(len(sightline.peaks))
            + ' lines are better than '
            + str(len(sightline.peaks) - 1)
        )
    else:
        print(
            str(len(sightline.peaks))
            + ' lines are NOT better than '
            + str(len(sightline.peaks) - 1)
        )
        sys.exit()

    # %%%%%%%%%%%%%%%%%%%%%%%%% Calculate likelihoods %%%%%%%%%%%%%%%%%%%%%%%%%
    likelihoods = []
    bary_likelihoods = []
    for i in range(len(sightlines)):
        resid = telluric_resids[i]
        bary_resid = bary_resids[i]
        sightline = sightlines[i]

        log_likelihood = np.log(1 / (sightline.sigma * np.sqrt(2 * np.pi))) + (
            -((resid) ** 2 / (2 * sightline.sigma ** 2))
        )

        bary_log_likelihood = np.log(1 / (sightline.sigma * np.sqrt(2 * np.pi))) + (
            -((bary_resid) ** 2 / (2 * sightline.sigma ** 2))
        )

        likelihoods.append(log_likelihood)
        bary_likelihoods.append(bary_log_likelihood)

    # %%%%%%%%%%%%%%%%%%%%%%% Co-add residuals and plot %%%%%%%%%%%%%%%%%%%%%%%
    coadd = np.ones_like(sightlines[0].grid)
    bary_coadd = np.ones_like(bary_grid)

    for i in range(len(sightlines)):
        likelihood = likelihoods[i]
        bary_likelihood = bary_likelihoods[i]

        coadd *= likelihood
        bary_coadd *= bary_likelihood

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare co-adds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if np.min(coadd) < np.min(bary_coadd):
        next_line = "Telluric"
    else:
        next_line = "Barycentric"

    name_next = "line" + str(iteration)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add new line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if next_line == "Telluric":
        idx = np.argmin(coadd)

        for sightline in sightlines:
            assert np.min(sightline.interp_flux[idx]) == sightline.interp_flux[idx]

        for i in range(len(sightlines)):
            sightline = sightlines[i]

            # get position and tau_0 of new line
            lam_0 = sightline.grid[idx]
            if len(sightline.peaks) > 3:
                tau_0 = 0.02
            else:
                tau_0 = 0.2
            line = {"lam_0": lam_0, "d": 0.01, "tau_0": tau_0}
            sightline.add_line(name=name_next, source=next_line, pars=line)

    elif next_line == "Barycentric":

        cloud_name = 'Cloud_' + str(sightline.num_sources)

        idx = np.argmin(bary_coadd)

        for i in range(len(sightlines)):
            sightline = sightlines[i]

            # get position and tau_0 of new line
            bary_lam_0 = bary_grid[idx]
            lam_0 = (
                bary_lam_0 - (sightline.v_bary / cst.c.to("km/s").value) * bary_lam_0
            )

            if len(sightline.peaks) > 3:
                tau_0 = 0.02
            else:
                tau_0 = 0.2
            line = {"lam_0": lam_0, "d": 0.001, "tau_0": tau_0}
            sightline.add_line(name=name_next, source=cloud_name, pars=line)

    else:
        print("Something bad happened...")
        sys.exit()

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot:

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Fit data
        fig, axs = plt.subplots(5, 2, sharex=True)
        for i in range(len(sightlines)):
            sightline = sightlines[i]
            out = outs[i]
            bary_out_interp = bary_outs[i]
            telluric_resid = telluric_resids[i]
            bary_resid_interp = bary_resids[i]

            for peak in sightline.peaks:
                if 'Telluric' in peak.name:
                    c = 'g'
                    label = 'Telluric lines'
                else:
                    c = 'r'
                    label = 'Barycentric lines'

                if sightline.most_recent in peak.name:
                    c = 'm'
                    label = 'Next line - not fit yet'

                axs[i, 0].vlines(
                    peak.value,
                    ymin=np.min(sightline.interp_flux),
                    ymax=np.max(sightline.interp_flux),
                    color=c,
                    label=label
                )

            axs[i, 0].plot(sightline.grid, sightline.interp_flux)
            axs[i, 0].plot(sightline.grid, out)
            axs[i, 0].plot(sightline.grid, telluric_resid)
            axs[i, 0].plot(sightline.grid, sightline.result.residual)

            axs[i, 1].plot(sightline.grid, sightline.interp_bary_flux)
            axs[i, 1].plot(bary_grid, bary_out_interp)
            axs[i, 1].plot(bary_grid, bary_resid_interp)
        plt.show()

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Likelihoods
        fig, axs = plt.subplots(1, 2)
        for i in range(len(likelihoods)):
            likelihood = likelihoods[i]
            axs[0].plot(sightline.grid, likelihood)
            axs[0].set_xlabel("Geocentric")

            bary_likelihood = bary_likelihoods[i]
            axs[1].plot(bary_grid, bary_likelihood)
            axs[1].set_xlabel("Barycentric")

        plt.suptitle("Likelihoods in each sightline")
        plt.show()

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Residuals
        fig, axs = plt.subplots(1, 2)

        axs[0].plot(sightlines[0].grid, coadd, label="Geocentric")
        axs[1].plot(bary_grid, bary_coadd, label="Barycentric")
        plt.suptitle("coadded likelihoods")
        plt.legend()
        plt.show()


def stop(sightlines):


    # get 'answer' in each sightline
    for i in range(len(sightlines)):
        sightline = sightlines[i]
        num_lines = len(sightline.peaks)

        # compute parameter ranges for each line and multiply
        range_terms = []
        for i in range(num_lines):

            lam_name = sightline.peaks[i].name
            b_name = lam_name.replace('lam_0', 'b')
            d_name = lam_name.replace('lam_0', 'd')
            tau_name = lam_name.replace('lam_0', 'tau_0')
            # print(lam_name, b_name, d_name, tau_name)

            lam_range = sightline.result.params[lam_name].max - sightline.result.params[lam_name].min
            b_range = sightline.result.params[b_name].max - sightline.result.params[b_name].min
            d_range = sightline.result.params[d_name].max - sightline.result.params[d_name].min
            tau_range = sightline.result.params[tau_name].max - sightline.result.params[tau_name].min

            range_term = lam_range * b_range * d_range * tau_range
            range_terms.append(range_term)

        # compute range for continuum parameters
        c_ranges = []
        for name in sightline.model_pars:
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

        a = np.math.factorial(num_lines) / denominator_term
        b = (4 * np.pi) ** ((2 * num_lines + sightline.n_anchors) / 2) / (np.sqrt(det))
        c = np.exp(sightline.result.chisqr / 2)
        answer = a * b * c

        try:
            sightline.answers.append(answer)
        except AttributeError:
            sightline.answers = [answer]

        print()
        print(sightline.datetime.date())
        print('a: ', a)
        print('b: ', b)
        print('a * b: ', a * b)
        print('c: ', c)
        print('answer: ', answer)
        print(sightline.answers)


for i in range(12):
    # print(sightlines[0].model_pars)

    loop(sightlines, i, method='least_squares', debug=True, plot=True)
