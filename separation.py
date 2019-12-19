from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy import constants as cst
import pandas as pd

from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.fit.models.create_model import createCont
from edibles.edibles.fit.models.models import Sightline
from edibles.edibles.fit.fit import fit

from sourceseparation.wavelength_corr import correctWavelength


def normalize(data, star_name, year):

    prominence = (np.max(data[1]) - np.min(data[1])) * 0.1
    peaks, _ = find_peaks(-data[1], prominence=prominence)
    peak_wavelengths = data[0].iloc[peaks]

    cont = createCont(data, n_points=4)

    sightline = Sightline(star_name=star_name, cont=cont)
    for i in range(len(peak_wavelengths)):
        source_name = "Source " + str(i + 1)
        line_name = "line" + str(i + 1)

        sightline.addSource(source_name=source_name, b=1.0, d=0.04)
        sightline.addLine(name=line_name, lam_0=peak_wavelengths.iloc[i], tau_0=0.5)

    fit_m, params = fit(star_name, data, sightline.model, breakdown=True, silent=True)

    f_norm = data[1] / cont(data[0])
    norm_data = [data[0], f_norm]

    return norm_data


def separate(observations, stop, xmin, xmax, flats, shifts):
    """

    :param observations: Set of n spectrum objects
    :type observations: list
    :param stop: Stopping condition
    :type stop: int (so far)
    :param xmin: Minimum wavelength value to fit
    :type xmin: float
    :param xmax: Maximum wavelength value to fit
    :type xmax: float
    :param flats: Wavelength range of flat section - to get sigma
    :type flats: tuple
    :param shifts: If true, calculate wavelength shifts for telluric lines
    :type shifts: bool


    :return: sightline objects with models fit to data
    :rtype: list
    
    """

    sigmas = []
    sightlines = []

    flat_xmin, flat_xmax = flats

    if shifts == True:
        shifts = correctWavelength(
            observations, xmin, xmax, silent_fit=True, silent_plot=True
        )
        print("Done calculating shifts")

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # INITIALIZATION

    datas_telluric = []
    datas_bary = []

    for i in range(len(observations)):

        obs = observations[i]

        # create wavelength data for each observation
        df_subset = obs.getSpectrum(xmin=xmin, xmax=xmax)
        data_telluric = (df_subset["wave"], df_subset["flux"])

        if shifts == True:
            # shift data
            data_telluric = (np.add(data_telluric[0], shifts[i]), data_telluric[1])

        # could streamline this a bit...
        data_telluric = normalize(data_telluric, obs.target, obs.date[0:4])

        data_bary = (df_subset["bary_wave"], data_telluric[1])
        datas_telluric.append(data_telluric)
        datas_bary.append(data_bary)

    # create sightline object
    for i in range(len(observations)):


        obs = observations[i]
        data_telluric = datas_telluric[i]
        data_bary = datas_bary[i]

        df_subset_flat = obs.getSpectrum(flat_xmin, flat_xmax)
        flat_data = (df_subset_flat["wave"], df_subset_flat["flux"])
        sigma_real = np.std(flat_data[1])
        mean = np.mean(flat_data[1])
        sigma = sigma_real / mean
        sigmas.append(sigma)


        cont = createCont(data_telluric, n_points=3)
        sightline = Sightline(star_name=obs.target, cont=cont)
        sightline.addSource(source_name="Telluric", b=2, d=0.01)
        sightline.addSource(source_name="Interstellar", b=1, d=0.001)
        sightlines.append(sightline)


    plt.figure()
    for data_telluric in datas_telluric:
        plt.plot(data_telluric[0], data_telluric[1])
    plt.title("Overplot of all telluric reference frames")
    # plt.show()

    print("Finished initializing")

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Fitting loop
    done = False
    counter = 0
    tau = 0.5


    COLNAMES1 = ["wave1", "flux1", "model1", "likelihood1"]
    COLNAMES2 = ["wave2", "flux2", "model2", "likelihood2"]
    COLNAMES3 = ["wave3", "flux3", "model3", "likelihood3"]
    COLNAMES4 = ["wave4", "flux4", "model4", "likelihood4"]
    COLNAMES5 = ["wave5", "flux5", "model5", "likelihood5"]

    colnames_list = [COLNAMES1, COLNAMES2, COLNAMES3, COLNAMES4, COLNAMES5]

    while not done:

        df_night1 = pd.DataFrame(columns=COLNAMES1)
        df_night2 = pd.DataFrame(columns=COLNAMES2)
        df_night3 = pd.DataFrame(columns=COLNAMES3)
        df_night4 = pd.DataFrame(columns=COLNAMES4)
        df_night5 = pd.DataFrame(columns=COLNAMES5)
        nights = [df_night1, df_night2, df_night3, df_night4, df_night5]

        df_night_bary1 = pd.DataFrame(columns=COLNAMES1)
        df_night_bary2 = pd.DataFrame(columns=COLNAMES2)
        df_night_bary3 = pd.DataFrame(columns=COLNAMES3)
        df_night_bary4 = pd.DataFrame(columns=COLNAMES4)
        df_night_bary5 = pd.DataFrame(columns=COLNAMES5)
        nights_bary = [
            df_night_bary1,
            df_night_bary2,
            df_night_bary3,
            df_night_bary4,
            df_night_bary5,
        ]

        for i in range(len(observations)):
            obs = observations[i]
            sigma = sigmas[i]
            sightline = sightlines[i]
            data_telluric = datas_telluric[i]
            data_bary = datas_bary[i]
            df_night = nights[i]
            df_night_bary = nights_bary[i]
            colnames = colnames_list[i]

            # Wave
            df_night[colnames[0]] = data_telluric[0]
            df_night[colnames[1]] = data_telluric[1]
            # Flux
            df_night_bary[colnames[0]] = data_bary[0]
            df_night_bary[colnames[1]] = data_bary[1]
            # Model
            fit_model = fit(
                obs.target,
                (df_night[colnames[0]], df_night[colnames[1]]),
                sightline.model,
                silent=True,
            )
            model_tell = fit_model(df_night[colnames[0]])
            df_night[colnames[2]] = model_tell
            df_night_bary[colnames[2]] = model_tell

        lefts = [
            df_night_bary1.iloc[0, 0],
            df_night_bary2.iloc[0, 0],
            df_night_bary3.iloc[0, 0],
            df_night_bary4.iloc[0, 0],
            df_night_bary5.iloc[0, 0],
        ]
        left = np.max(lefts)
        rights = [
            df_night_bary1.iloc[-1, 0],
            df_night_bary2.iloc[-1, 0],
            df_night_bary3.iloc[-1, 0],
            df_night_bary4.iloc[-1, 0],
            df_night_bary5.iloc[-1, 0],
        ]
        right = np.min(rights)

        df_night1 = df_night1[df_night1["wave1"].between(left, right)]
        df_night2 = df_night2[df_night2["wave2"].between(left, right)]
        df_night3 = df_night3[df_night3["wave3"].between(left, right)]
        df_night4 = df_night4[df_night4["wave4"].between(left, right)]
        df_night5 = df_night5[df_night5["wave5"].between(left, right)]
        nights = [df_night1, df_night2, df_night3, df_night4, df_night5]

        df_night_bary1 = df_night_bary1[df_night_bary1["wave1"].between(left, right)]
        df_night_bary2 = df_night_bary2[df_night_bary2["wave2"].between(left, right)]
        df_night_bary3 = df_night_bary3[df_night_bary3["wave3"].between(left, right)]
        df_night_bary4 = df_night_bary4[df_night_bary4["wave4"].between(left, right)]
        df_night_bary5 = df_night_bary5[df_night_bary5["wave5"].between(left, right)]
        nights_bary = [
            df_night_bary1,
            df_night_bary2,
            df_night_bary3,
            df_night_bary4,
            df_night_bary5,
        ]

        for i in range(len(nights)):
            sigma = sigmas[i]
            df_night = nights[i]
            df_night_bary = nights_bary[i]
            colnames = colnames_list[i]

            # Likelihood
            resid = df_night[colnames[1]] - df_night[colnames[2]]

            log_likelihood_telluric = np.log(1 / (sigma * np.sqrt(2 * np.pi))) + (
                -((resid) ** 2 / (2 * sigma ** 2))
            )
            df_night[colnames[3]] = log_likelihood_telluric

            # like_t = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(((resid) ** 2 / (2 * sigma**2)))
            # df_night[colnames[3]] = like_t

            resid_bary = df_night_bary[colnames[1]] - df_night_bary[colnames[2]]

            log_likelihood_bary = np.log(1 / (sigma * np.sqrt(2 * np.pi))) + (
                -((resid_bary) ** 2 / (2 * sigma ** 2))
            )
            df_night_bary[colnames[3]] = log_likelihood_bary

            # like_b = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(((resid_bary) ** 2 / (2 * sigma**2)))
            # df_night_bary[colnames[3]] = like_b

        for i in range(len(nights)):

            nights[i] = nights[i].dropna().reset_index()
            nights_bary[i] = nights_bary[i].dropna().reset_index()

        df_tell = pd.concat(
            [nights[0], nights[1], nights[2], nights[3], nights[4]], axis=1, sort=False
        )
        df_bary = pd.concat(
            [
                nights_bary[0],
                nights_bary[1],
                nights_bary[2],
                nights_bary[3],
                nights_bary[4],
            ],
            axis=1,
            sort=False,
        )

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # PLOTTING
        fig, axes = plt.subplots(nrows=5, ncols=2, sharex="col", sharey=False)
        title = "Iteration " + str(counter + 1) + ", data, model and residual"
        fig.suptitle(title)
        i = 0
        for row in axes:
            j = 0
            for col in row:
                if j == 0:
                    # data

                    colnames = colnames_list[i]

                    col.plot(df_tell[colnames[0]], df_tell[colnames[1]])
                    # model
                    col.plot(df_tell[colnames[0]], df_tell[colnames[2]])
                    # residual
                    resid = df_tell[colnames[1]] - df_tell[colnames[2]]
                    col.plot(df_tell[colnames[0]], resid)
                    # one_lam_val = sightlines[0].lines['all'][-1].lam_0.val
                    # col.axvline(x=one_lam_val, c='r', label=str(one_lam_val))

                    if i == 4:
                        col.set_xlabel("Telluric frame")
                    ylabel = "Obs #" + str(i)
                    col.set_ylabel(ylabel)
                elif j == 1:
                    # data
                    col.plot(df_bary[colnames[0]], df_bary[colnames[1]])
                    # model
                    col.plot(df_bary[colnames[0]], df_bary[colnames[2]])
                    # residual
                    resid = df_bary[colnames[1]] - df_bary[colnames[2]]
                    col.plot(df_bary[colnames[0]], resid)

                    # col.axvline(x=7664.55, c='r', label=str(7664.55))

                    if i == 4:
                        col.set_xlabel("Barycentric frame")
                j += 1
            i += 1
        # plt.show()

        fig, axes = plt.subplots(nrows=6, ncols=2, sharex="col", sharey=False)
        title = "Iteration " + str(counter + 1) + ", likelihoods"
        fig.suptitle(title)
        i = 0
        for i in range(5):
            row = axes[i]
            j = 0
            for col in row:
                if j == 0:
                    colnames = colnames_list[i]
                    col.plot(df_tell[colnames[0]], df_tell[colnames[3]])
                    axes[5, 0].plot(df_tell[colnames[0]], df_tell[colnames[3]])

                    if i == 4:
                        col.set_xlabel("Telluric frame")
                    ylabel = "Obs #" + str(i)
                    col.set_ylabel(ylabel)
                elif j == 1:
                    colnames = colnames_list[i]
                    col.plot(df_bary[colnames[0]], df_bary[colnames[3]])
                    axes[5, 1].plot(df_bary[colnames[0]], df_bary[colnames[3]])

                    if i == 4:
                        col.set_xlabel("Barycentric frame")
                j += 1
            i += 1

        # plt.show()

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # calculate combined likelihood

        # TELLURIC
        combined = (
            df_tell["likelihood1"]
            * df_tell["likelihood2"]
            * df_tell["likelihood3"]
            * df_tell["likelihood4"]
            * df_tell["likelihood5"]
        )

        idx_tell = np.argmin(combined)
        guess_wave_tell = nights[0][colnames_list[0][0]][idx_tell]
        guess_flux_tell = combined[idx_tell]

        # BARYCENTRIC
        combined_bary = (
            df_bary["likelihood1"]
            * df_bary["likelihood2"]
            * df_bary["likelihood3"]
            * df_bary["likelihood4"]
            * df_bary["likelihood5"]
        )

        idx_bary = np.argmin(combined_bary)
        guess_wave_bary = nights_bary[0][colnames_list[0][0]][idx_bary]
        guess_flux_bary = combined_bary[idx_bary]

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # PLOTTING
        fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True)
        title = "Iteration " + str(counter + 1) + ", combined likelihoods"
        fig.suptitle(title)
        # TELLURIC
        axes[0].plot(df_tell["wave1"], combined, "k")
        axes[0].plot(guess_wave_tell, guess_flux_tell, "x")
        axes[0].set_title("Telluric")

        # BARYCENTRIC
        axes[1].plot(df_bary["wave1"], combined_bary, "k")
        axes[1].plot(guess_wave_bary, guess_flux_bary, "x")
        axes[1].set_title("Barycentric")

        # plt.show()


        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # choose ref frame to make peak

        # TELLURIC
        if guess_flux_tell < guess_flux_bary:
            telluric_ref_wave = guess_wave_tell

            # TODO: implement something for oxygen lines - linked wavelength for 2 peaks

            turn = "Telluric"
            info = "Iteration: " + str(counter) + "     Next line is: " + turn
            print(info)

            for i in range(len(observations)):
                obs = observations[i]
                sightline = sightlines[i]

                sightline.setSource(source_name="Telluric")
    
                print()
                print("Observation #", str(i))
                print("Lines so far:")
                line_names = [line.name for line in sightline.lines["all"]]
                line_vals = [str(line.lam_0.val) for line in sightline.lines["all"]]
                for i in range(len(line_names)):
                    print(line_names[i] + ":  " + line_vals[i])

                new_line_name = "tell_" + str(counter)
                print("New line: " + new_line_name + " -> " + str(telluric_ref_wave))
                sightline.addLine(
                    name=new_line_name, lam_0=telluric_ref_wave, tau_0=tau
                )

        # BARYCENTRIC
        elif guess_flux_bary < guess_flux_tell:
            bary_ref_wave = guess_wave_bary

            turn = "Interstellar"
            info = "Iteration: " + str(counter) + "     Next line is: " + turn
            print(info)

            for i in range(len(observations)):
                obs = observations[i]
                sightline = sightlines[i]


                telluric_ref_wave = bary_ref_wave / (
                    1 + (obs.v_bary / cst.c.to("km/s").value)
                )

                sightline.setSource(source_name="Interstellar")

                print()
                print("Observation #", str(i))
                print("Lines so far:")
                line_names = [line.name for line in sightline.lines["all"]]
                line_vals = [str(line.lam_0.val) for line in sightline.lines["all"]]
                for i in range(len(line_names)):
                    print(line_names[i] + ":  " + line_vals[i])

                new_line_name = "int_" + str(counter)
                print("New line: " + new_line_name + " -> " + str(telluric_ref_wave))
                sightline.addLine(
                    name=new_line_name, lam_0=telluric_ref_wave, tau_0=tau
                )
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        print()
        print()
        print()

        counter += 1
        if counter >= stop:
            done = True


        # update tau for next loop
        print(sightlines[0].lines["all"])
        newest_line = sightlines[0].lines["all"][-1]
        tau = newest_line.tau_0.val

    return sightlines

if __name__ == "__main__":
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # HD170740 KI test region

    file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"
    sp1 = EdiblesSpectrum(file1)
    sp2 = EdiblesSpectrum(file2)
    sp3 = EdiblesSpectrum(file3)
    sp4 = EdiblesSpectrum(file4)
    sp5 = EdiblesSpectrum(file5)

    observations = [sp1, sp2, sp3, sp4, sp5]

    # xmin = 7661.5
    # xmax = 7669.0
    xmin = 7658
    xmax = 7675

    flat_xmin = 7667.6
    flat_xmax = 7669.1
    flats = (flat_xmin, flat_xmax)

    sightlines = separate(
        observations, stop=10, xmin=xmin, xmax=xmax, flats=flats, shifts=True
    )

    print(sightlines)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # # HD170740 NaI test region

    # file1 = '/HD170740/BLUE_346/HD170740_w346_blue_20140916_O12.fits'
    # file2 = '/HD170740/BLUE_346/HD170740_w346_blue_20150424_O12.fits'
    # file3 = '/HD170740/BLUE_346/HD170740_w346_blue_20160505_O12.fits'
    # file4 = '/HD170740/BLUE_346/HD170740_w346_blue_20160612_O12.fits'
    # file5 = '/HD170740/BLUE_346/HD170740_w346_blue_20170701_O12.fits'
    # sp1 = EdiblesSpectrum(file1)
    # sp2 = EdiblesSpectrum(file2)
    # sp3 = EdiblesSpectrum(file3)
    # sp4 = EdiblesSpectrum(file4)
    # sp5 = EdiblesSpectrum(file5)

    # observations = [sp1, sp2, sp3, sp4, sp5]
    # xmin = 3301
    # xmax = 3304
    # flat_xmin = 3303
    # flat_xmax = 3303.5
    # flats = (flat_xmin, flat_xmax)

    # sightlines = separate(observations, stop=3, xmin=xmin, xmax=xmax, flats=flats shifts=False)

    # print(sightlines)

