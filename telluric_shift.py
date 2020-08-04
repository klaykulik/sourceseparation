import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

from edibles.sightline import Sightline

from sourceseparation.read_hitran import read_hitran, convert


def telluric_shift(observations, xmin, xmax, zoom_min, zoom_max, plot=False):
    '''Correct the spectra using telluric lines from the HITRAN dataset.

    xmin/xmax are the wavelength bounds of the data to fit telluric lines to.
    zoom_min/zoom_max are the bounds for the returned EdiblesSpectrum objects.


    '''

    assert zoom_min > xmin + 0.5
    assert zoom_max < xmax - 0.5

    # O_2 data from HITRAN
    pars_list = convert(read_hitran('telluric_lines_HITRAN.txt'))

    # Set tau_cutoff (0.4, 0.02, or 0.0)
    tau_cutoff = 0.4

    linelist = []
    for pars in pars_list:
        if (xmin < pars['lam_0']) and (pars['lam_0'] < xmax):
            if pars['tau_0'] > tau_cutoff:
                linelist.append(pars)
    linelist.reverse()


    sightlines = []
    for sp in observations:
        sp.getSpectrum(xmin, xmax)

        sightline = Sightline(sp, n_anchors=4)
        sightlines.append(sightline)

    resids = []
    outs = []
    errs_list = []
    lams_list = []
    for sightline in sightlines:
        print('Fitting sightline...')
        sightline.add_source(name='O2', similar={'b': 0.6})

        for i in range(len(linelist)):
            line = linelist[i]
            if line['tau_0'] > tau_cutoff:
                name = 'line' + str(i)
                sightline.add_line(name=name, source='O2', pars=line)
                par_name_d = 'O2_' + name + '_d'
                sightline.model_pars[par_name_d].set(value=0.05)

        sightline.fit(report=False, plot=False)

        out = sightline.model.eval(data=sightline.flux,
                                   params=sightline.result.params,
                                   x=sightline.wave)
        resid = sightline.flux - out
        resids.append(resid)
        outs.append(out)

        errs = []
        lams = []
        for name in sightline.model.param_names:
            if name[-5:] == 'lam_0':
                errs.append(sightline.result.params[name].stderr)
                lams.append(sightline.result.params[name].value)
        errs_list.append(errs)
        lams_list.append(lams)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # CORRECT WAVELENGTH SHIFT FROM FIT
    reals = [linelist[i]['lam_0'] for i in range(len(linelist))]
    shift_resids = [np.subtract(reals, lams) for lams in lams_list]

    cs = ['C0', 'C1', 'C2', 'C3', 'C4']

    shifts = []
    for i in range(len(observations)):
        sp = observations[i]
        lams = lams_list[i]
        errs = errs_list[i]
        shift_resid = shift_resids[i]

        # weights = np.mean(shift_resids[i]) - np.abs(np.mean(shift_resids[i]) - shift_resids[i])
        # weights = 1 / shift_resids[i]
        spl = UnivariateSpline(reals, shift_resids[i])
        # spl.set_smoothing_factor(5)

        if plot:
            plt.errorbar(reals, shift_resid, yerr=errs, fmt='o', c=cs[i],
                         label=sp.datetime.date())
            plt.plot(reals, spl(reals), linestyle='dashed')
            plt.plot(sp.wave, spl(sp.wave), c=cs[i])

        shift = spl(sp.wave)
        shifts.append(shift)

    if plot:
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Residual fit difference (Angstroms)')
        plt.title("HD170740")
        plt.legend()
        plt.show()

    for i in range(len(observations)):
        sp = observations[i]
        shift = shifts[i]
        sp.shift(shift, xmin=zoom_min, xmax=zoom_max)

    if plot:
        for sp in observations:
            plt.plot(sp.wave, sp.flux)
        plt.scatter(reals, np.zeros_like(reals), c='k', label='HITRAN data')
        plt.show()

    return observations


if __name__ == '__main__':

    import os

    from edibles.utils.edibles_spectrum import EdiblesSpectrum

    CURRENT_DIR = os.getcwd()

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    os.chdir(CURRENT_DIR)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # CORRECT TELLURIC SHIFT

    observations = telluric_shift(observations, xmin, xmax)








