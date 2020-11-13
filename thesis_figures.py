import matplotlib.pyplot as plt
import numpy as np

from edibles.utils.edibles_spectrum import EdiblesSpectrum


def reference_frames():

    file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
    file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

    xmin = 7661.5
    xmax = 7669

    # GEOCENTRIC
    sp1 = EdiblesSpectrum(file1)
    sp1.getSpectrum(xmin=xmin, xmax=xmax)
    sp1.flux = sp1.flux / np.max(sp1.flux)
    plt.plot(sp1.wave, sp1.flux, "k")
    plt.text(7661.6, 0.8, sp1.datetime.date(), fontsize=12)

    sp2 = EdiblesSpectrum(file2)
    sp2.getSpectrum(xmin=xmin, xmax=xmax)
    sp2.flux = sp2.flux / np.max(sp2.flux) + 1
    plt.plot(sp2.wave, sp2.flux, "k")
    plt.text(7661.6, 1.8, sp2.datetime.date(), fontsize=12)

    sp3 = EdiblesSpectrum(file3)
    sp3.getSpectrum(xmin=xmin, xmax=xmax)
    sp3.flux = sp3.flux / np.max(sp3.flux) + 2
    plt.plot(sp3.wave, sp3.flux, "k")
    plt.text(7661.6, 2.8, sp3.datetime.date(), fontsize=12)

    sp4 = EdiblesSpectrum(file4)
    sp4.getSpectrum(xmin=xmin, xmax=xmax)
    sp4.flux = sp4.flux / np.max(sp4.flux) + 3
    plt.plot(sp4.wave, sp4.flux, "k")
    plt.text(7661.6, 3.8, sp4.datetime.date(), fontsize=12)

    sp5 = EdiblesSpectrum(file5)
    sp5.getSpectrum(xmin=xmin, xmax=xmax)
    sp5.flux = sp5.flux / np.max(sp5.flux) + 4
    plt.plot(sp5.wave, sp5.flux, "k")
    plt.text(7661.6, 4.8, sp5.datetime.date(), fontsize=12)

    plt.title("HD 170740", fontsize=14)
    plt.xlabel(r"Wavelength ($\AA$)", fontsize=14)
    plt.ylabel("Normalized Flux + Offset", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

    # BARYCENTRIC
    sp1.bary_flux = sp1.bary_flux / np.max(sp1.bary_flux)
    plt.plot(sp1.bary_wave, sp1.bary_flux, "k")
    plt.text(7661.6, 0.8, sp1.datetime.date(), fontsize=12)

    sp2.bary_flux = sp2.bary_flux / np.max(sp2.bary_flux) + 1
    plt.plot(sp2.bary_wave, sp2.bary_flux, "k")
    plt.text(7661.6, 1.8, sp2.datetime.date(), fontsize=12)

    sp3.bary_flux = sp3.bary_flux / np.max(sp3.bary_flux) + 2
    plt.plot(sp3.bary_wave, sp3.bary_flux, "k")
    plt.text(7661.6, 2.8, sp3.datetime.date(), fontsize=12)

    sp4.bary_flux = sp4.bary_flux / np.max(sp4.bary_flux) + 3
    plt.plot(sp4.bary_wave, sp4.bary_flux, "k")
    plt.text(7661.6, 3.8, sp4.datetime.date(), fontsize=12)

    sp5.bary_flux = sp5.bary_flux / np.max(sp5.bary_flux) + 4
    plt.plot(sp5.bary_wave, sp5.bary_flux, "k")
    plt.text(7661.6, 4.8, sp5.datetime.date(), fontsize=12)

    plt.title("HD 170740", fontsize=14)
    plt.xlabel(r"Wavelength ($\AA$)", fontsize=14)
    plt.ylabel("Normalized Flux + Offset", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()


def line_examples():
    file = "/HD170740/RED_860/HD170740_w860_redl_20160613_O6.fits"

    xmin = 7055
    xmax = 7130

    # GEOCENTRIC
    sp1 = EdiblesSpectrum(file)
    sp1.getSpectrum(xmin=xmin, xmax=xmax)
    sp1.flux = sp1.flux / np.max(sp1.flux)
    plt.plot(sp1.wave, sp1.flux, "k")
    plt.title("HD 170740", fontsize=14)
    plt.xlabel(r"Wavelength ($\AA$)", fontsize=14)
    plt.ylabel("Normalized Flux", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()


def ica():

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import signal

    from sklearn.decomposition import FastICA

    # #############################################################################
    # Generate sample data
    np.random.seed(0)
    n_samples = 2000
    time = np.linspace(0, 8, n_samples)

    s1 = 2 * np.sin(2 * time)  # Signal 1 : sinusoidal signal
    s2 = np.sign(np.sin(3 * time))  # Signal 2 : square signal
    s3 = signal.sawtooth(2 * np.pi * time)  # Signal 3: saw tooth signal

    S = np.c_[s1, s2, s3]
    S += 0.2 * np.random.normal(size=S.shape)  # Add noise

    # S /= S.std(axis=0)  # Standardize data
    # Mix data
    A = np.array([[1, 1, 1], [0.5, 2, 1.0], [1.5, 1.0, 2.0]])  # Mixing matrix
    X = np.dot(S, A.T)  # Generate observations

    # Compute ICA
    ica = FastICA(n_components=3)
    S_ = ica.fit_transform(X)  # Reconstruct signals
    A_ = ica.mixing_  # Get estimated mixing matrix

    # We can `prove` that the ICA model applies by reverting the unmixing.
    assert np.allclose(X, np.dot(S_, A_.T) + ica.mean_)


    # #############################################################################
    # Plot results

    plt.figure()

    models = [X, S, S_]
    names = ['Observations (mixed signal)',
             'True Sources',
             'ICA recovered signals',
             'PCA recovered signals']
    colors = ['C0', 'C1', 'C2']

    for ii, (model, name) in enumerate(zip(models, names), 1):
        plt.subplot(3, 1, ii)
        plt.title(name)
        for sig, color in zip(model.T, colors):
            plt.plot(sig, color=color)

    plt.tight_layout()
    plt.show()


def bayes():

    import numpy as np
    from matplotlib import pyplot as plt

    import scipy.stats as stats



    A = 0.5




    dist = stats.beta
    n_trials = [0, 1, 2, 3, 4, 100]
    data = stats.bernoulli.rvs(A, size=n_trials[-1])
    x = np.linspace(0, 1, 100)
    prior_list = []
    heads_list = []
    # For the already prepared, I'm using Binomial's conj. prior.
    # for k, N in enumerate(n_trials):
    #     heads = data[:N].sum()
    #     prior = dist.pdf(x, 1 + heads, 1 + N - heads)

    #     prior_list.append(prior)
    #     heads_list.append(heads)


    # for i in range(len(prior_list)):


        # if i > 0:
        #     # prior
        #     plt.plot(x, prior_list[i - 1], label="Prior")
        #     plt.fill_between(x, 0, prior_list[i - 1], color="C0", alpha=0.25)
        #     plt.vlines(0.5, 0, 4, color="k", linestyles="--", lw=1)
        #     plt.xlabel("$p$, probability of heads", fontsize=14)
        #     plt.xticks(fontsize=12)
        #     plt.yticks(fontsize=12)
        #     plt.autoscale(tight=True)
        #     plt.show()

        #     # likelihood
        #     plt.plot(x, prior_list[i], label="Likelihood")
        #     plt.fill_between(x, 0, prior_list[i], color="C1", alpha=0.25)
        #     plt.vlines(0.5, 0, 4, color="k", linestyles="--", lw=1)
        #     plt.xlabel("$p$, probability of heads", fontsize=14)
        #     plt.xticks(fontsize=12)
        #     plt.yticks(fontsize=12)
        #     plt.autoscale(tight=True)

        #     # posterior
        #     plt.plot(x, prior_list[i] * prior_list[i - 1], linestyle='--', label="Posterior")
        #     plt.fill_between(x, 0, prior_list[i] * prior_list[i - 1], color="C2", alpha=0.25)
        #     plt.vlines(0.5, 0, 4, color="k", linestyles="--", lw=1)
        #     plt.xlabel("$p$, probability of heads", fontsize=14)
        #     plt.xticks(fontsize=12)
        #     plt.yticks(fontsize=12)
        #     plt.autoscale(tight=True)


        #     plt.legend(fontsize=14)
        #     plt.show()


    # For the already prepared, I'm using Binomial's conj. prior.
    for k, N in enumerate(n_trials):
        sx = plt.subplot(len(n_trials)/2, 2, k+1)
        plt.xlabel("Probability of heads") \
            if (k == 4 or k == 5) else None
        plt.setp(sx.get_yticklabels(), visible=False)
        heads = data[:N].sum()
        y = dist.pdf(x, 1 + heads, 1 + N - heads)
        plt.plot(x, y, label="observe %d tosses,\n %d heads" % (N, heads))
        plt.fill_between(x, 0, y, color="#348ABD", alpha=0.4)
        plt.vlines(0.5, 0, 4, color="k", linestyles="--", lw=1)

        leg = plt.legend()
        leg.get_frame().set_alpha(0.4)
        # plt.autoscale(tight=True)
    plt.tight_layout()
    plt.show()



def coadding():


    file1 = "/HD170740/BLUE_437/HD170740_w437_blue_20140915_O15.fits"
    file2 = "/HD170740/BLUE_437/HD170740_w437_blue_20140916_O15.fits"
    file3 = "/HD170740/BLUE_437/HD170740_w437_blue_20150626_O15.fits"
    file4 = "/HD170740/BLUE_437/HD170740_w437_blue_20160613_O15.fits"
    file5 = "/HD170740/BLUE_437/HD170740_w437_blue_20170705_O15.fits"

    xmin = 4220
    xmax = 4227.5

    sp1 = EdiblesSpectrum(file1)
    sp2 = EdiblesSpectrum(file2)
    sp3 = EdiblesSpectrum(file3)
    sp4 = EdiblesSpectrum(file4)
    sp5 = EdiblesSpectrum(file5)

    observations = [sp1, sp2, sp3, sp4, sp5]

    fig, axs = plt.subplots(1, 2)

    shift = 0
    for sp in observations:
        shift += 0.01
        sp.getSpectrum(xmin=xmin, xmax=xmax)
        sp.interp_flux = sp.interp_flux / np.max(sp.interp_flux)
        axs[0].plot(sp.grid, sp.interp_flux + shift)


    superspectrum = np.ones_like(observations[0].interp_flux)
    for sp in observations:
        superspectrum = superspectrum * sp.interp_flux

    axs[0].plot(observations[0].grid, superspectrum, 'k')
    axs[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=14)
    axs[0].set_title("Geocentric Reference Frame", fontsize=14)

    axs[0].tick_params(axis='both', labelsize=12)
    # axs[0].set_yticklabels(fontsize=12)



    shift = 0
    for sp in observations:
        shift += 0.01
        sp.getSpectrum(xmin=xmin, xmax=xmax)
        sp.interp_bary_flux = sp.interp_bary_flux / np.max(sp.interp_bary_flux)
        axs[1].plot(sp.grid, sp.interp_bary_flux + shift)


    superspectrum = np.ones_like(observations[0].interp_bary_flux)
    for sp in observations:
        superspectrum = superspectrum * sp.interp_bary_flux


    axs[1].plot(observations[0].grid, superspectrum, 'k')
    axs[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=14)
    axs[1].set_title("Barycentric Reference Frame", fontsize=14)


    axs[1].tick_params(axis='both', labelsize=12)
    # axs[0].set_yticklabels(fontsize=12)



    plt.show()



if __name__ == "__main__":

    # reference_frames()

    # line_examples()

    # ica()


    # bayes()

    coadding()
