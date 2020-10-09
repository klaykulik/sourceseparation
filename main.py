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

for sp in observations:
    # plt.plot(sp.wave, sp.flux)
    plt.plot(sp.grid, sp.interp_flux)
plt.show()

# Ensure all same length
length = [len(sp.interp_flux) for sp in observations]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE NEW MODELS

sightlines = []
for i in range(len(observations)):
    sp = observations[i]
    sightline = Sightline(sp, n_anchors=4)
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
flat_xmin = 7667
flat_xmax = 7668

for sightline in sightlines:
    sigma_idx = np.where(
        np.logical_and(sightline.grid > flat_xmin, sightline.grid < flat_xmax)
    )
    sightline.sigma = np.std(sightline.grid[sigma_idx])
    # print('Sigma: ', sightline.sigma)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LOOP


def loop(sightlines, iteration, debug=False, plot=False):

    # %%%%%%%%%%%%%%%% Fit geocentric model and get residuals %%%%%%%%%%%%%%%%%

    telluric_resids = []
    outs = []
    for i in range(len(sightlines)):
        sightline = sightlines[i]

        print("Fitting sightline...")
        sightline.fit(
            data=sightline.interp_flux, x=sightline.grid, report=debug, plot=debug
        )
        out = sightline.model.eval(
            data=sightline.interp_flux, params=sightline.result.params, x=sightline.grid
        )
        outs.append(out)

        telluric_resid = sightline.interp_flux - out
        telluric_resids.append(telluric_resid)

        print(sightline.result.covar)

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

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot:
        fig, axs = plt.subplots(5, 2, sharex=True)
        for i in range(len(sightlines)):
            sightline = sightlines[i]
            out = outs[i]
            bary_out_interp = bary_outs[i]
            telluric_resid = telluric_resids[i]
            bary_resid_interp = bary_resids[i]

            axs[i, 0].plot(sightline.grid, sightline.interp_flux)
            axs[i, 0].plot(sightline.grid, out)
            axs[i, 0].plot(sightline.grid, telluric_resid)

            axs[i, 1].plot(sightline.grid, sightline.interp_bary_flux)
            axs[i, 1].plot(bary_grid, bary_out_interp)
            axs[i, 1].plot(bary_grid, bary_resid_interp)
        plt.show()

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

    # %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot likelihoods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot:
        for likelihood in likelihoods:
            plt.plot(sightline.grid, likelihood)
        plt.title("Geocentric Likelihoods")
        plt.show()
        for bary_likelihood in bary_likelihoods:
            plt.plot(bary_grid, bary_likelihood)
        plt.title("Barycentric Likelihoods")
        plt.show()

    # %%%%%%%%%%%%%%%%%%%%%%% Co-add residuals and plot %%%%%%%%%%%%%%%%%%%%%%%
    coadd = np.ones_like(sightlines[0].grid)
    bary_coadd = np.ones_like(bary_grid)

    for i in range(len(sightlines)):
        likelihood = likelihoods[i]
        bary_likelihood = bary_likelihoods[i]

        coadd *= likelihood
        bary_coadd *= bary_likelihood

    if plot:
        fig, axs = plt.subplots(1, 2)

        axs[0].plot(sightlines[0].grid, coadd, label="Geocentric")
        axs[1].plot(bary_grid, bary_coadd, label="Barycentric")
        plt.suptitle("coadded likelihoods")
        plt.legend()
        plt.show()

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare co-adds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if np.min(coadd) < np.min(bary_coadd):
        next_line = "Telluric"
    else:
        next_line = "Barycentric"

    name = "line" + str(iteration)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add new line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if next_line == "Telluric":
        idx = np.argmin(coadd)

        for sightline in sightlines:
            assert np.min(sightline.interp_flux[idx]) == sightline.interp_flux[idx]

        for i in range(len(sightlines)):
            sightline = sightlines[i]

            # get position and tau_0 of new line
            lam_0 = sightline.grid[idx]
            tau_0 = 0.2
            line = {"lam_0": lam_0, "d": 0.01, "tau_0": tau_0}
            sightline.add_line(name=name, source=next_line, pars=line)

    elif next_line == "Barycentric":
        idx = np.argmin(bary_coadd)

        for i in range(len(sightlines)):
            sightline = sightlines[i]

            # get position and tau_0 of new line
            bary_lam_0 = bary_grid[idx]
            lam_0 = (
                bary_lam_0 - (sightline.v_bary / cst.c.to("km/s").value) * bary_lam_0
            )

            tau_0 = 0.2

            line = {"lam_0": lam_0, "d": 0.01, "tau_0": tau_0}
            sightline.add_line(name=name, source=next_line, pars=line)

    else:
        print("Something bad happened...")
        sys.exit()


answers = []
for i in range(7):
    loop(sightlines, i, debug=False, plot=True)

    # the equation:
    num_lines = len(sightlines[0].peaks)
    lam_range = 1
    b_range = 1
    d_range = 1
    tau_range = 1
    c_range = 1
    chi_sq_min = 1
    n_anchors = sightline.n_anchors
    det = np.linalg.det(sightline.result.covar)

    answer = (
        np.math.factorial(num_lines)
        / (
            (lam_range * b_range * d_range * tau_range) ** num_lines
            * c_range ** n_anchors
        )
        * (4 * np.pi) ** ((2 * num_lines + n_anchors) / 2)
        / (np.sqrt(det))
        * np.exp(chi_sq_min / 2)
    )


    print(answer)

    answers.append(answer)

    if len(answers) > 1:
        if answer / answers[-1] > 1:
            print(str(num_lines) + ' lines are better than ' + str(num_lines - 1))

        else:
            print('Not better!')
            sys.exit()
