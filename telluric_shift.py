import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

from edibles.sightline import Sightline
from edibles.utils.edibles_spectrum import EdiblesSpectrum

from sourceseparation.read_hitran import read_hitran, convert


def telluric_shift(observations, xmin, xmax, zoom_xmin, zoom_xmax, molecule='O2', plot=False):
    '''Correct the spectra using telluric lines from the HITRAN dataset.

    Args:
        observations (list): A lsit of EdiblesSpectrum objects
        xmin (float): Minimum wavelength to fit the telluric spectrum,
            There should be lines in this range from HITRAN
        xmax (float): Maximum wavelength to fit the telluric spectrum,
            There should be lines in this range from HITRAN
        zoom_xmin (float): Minimum wavelength of the returned spectrum,
            must be greater than xmin and xmin+shift
        zoom_xmax (float): Maximum wavelength of the returned spectrum,
            must be less than xmax and xmax+shift
        molecule (str): The molecule to get HITRAN data from
        plot (bool): If true, plots the fitting of each sightline,
            the peak differences from the HITRAN data (with spline), and the shifted data.

    Returns:
        list: The initial list of EdiblesSpectrum objects, each with updated (shifted)
            wave, bary_wave, flux, bary_flux, grid, interp_flux, and interp_bary_flux attributes

    '''

    assert zoom_xmin > xmin + 0.5
    assert zoom_xmax < xmax - 0.5


    # telluric data from HITRAN
    telluric_files_dict = {
        'O2': 'telluric_lines_HITRAN.txt',
        'H2O': 'telluric_H2O_lines_HITRAN.txt'
    }
    telluric_filename = telluric_files_dict[molecule]

    pars_list = convert(read_hitran(telluric_filename, molecule))

    # Set tau_cutoff (0.4, 0.02, or 0.0)
    tau_cutoff = 0.4

    linelist = []
    for pars in pars_list:
        if (xmin < pars['lam_0']) and (pars['lam_0'] < xmax):
            if pars['tau_0'] > tau_cutoff:
                if not (7663 < pars['lam_0'] < 7668):
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
        print('Fitting telluric lines and shifting wavelengths...')

        for i in range(len(linelist)):
            line = linelist[i]
            if line['tau_0'] > tau_cutoff:
                name = 'line' + str(i)
                line['d'] = 0.05
                sightline.add_line(name=name, source='Telluric', pars=line)
                # par_name_d = 'Telluric_' + name + '_d'
                # sightline.all_pars[par_name_d].set(value=0.05)

        sightline.fit(report=False, plot=plot)

        out = sightline.complete_model.eval(data=sightline.flux,
                                            params=sightline.result.params,
                                            x=sightline.wave)
        resid = sightline.flux - out
        resids.append(resid)
        outs.append(out)

        errs = []
        lams = []
        for name in sightline.complete_model.param_names:
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

        spl = UnivariateSpline(reals, shift_resid)

        if plot:
            plt.errorbar(reals, shift_resid, yerr=errs, fmt='o', c=cs[i],
                         label=sp.datetime.date())
            # plt.plot(reals, spl(reals), linestyle='dashed')
            plt.plot(sp.wave, spl(sp.wave), c=cs[i])

        shift = spl(sp.wave)
        shifts.append(shift)

    if plot:
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel(r'Difference between fit peaks and HITRAN peaks ($\AA$)')
        plt.title("HD170740")
        plt.legend()
        plt.show()

    zoom_reals = []
    for real in reals:
        if (zoom_xmin < real) and (real < zoom_xmax):
            zoom_reals.append(real)

    for i in range(len(observations)):
        sp = observations[i]
        shift = shifts[i]
        sp.shift(shift, zoom_xmin=zoom_xmin, zoom_xmax=zoom_xmax)

    if plot:
        for sp in observations:
            plt.plot(sp.wave, sp.flux / np.max(sp.flux), label=sp.datetime.date())
        plt.scatter(zoom_reals, np.zeros_like(zoom_reals), c='k', label='HITRAN data')
        plt.legend()
        plt.xlabel(r'Wavelength ($\AA$)', fontsize=14)
        plt.ylabel('Normalized Flux', fontsize=14)
        plt.tick_params(axis='both', labelsize=12)
        plt.show()

    return observations


if __name__ == '__main__':

    example = 1

    if example == 1:
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

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # CORRECT TELLURIC SHIFT

        observations = telluric_shift(observations, xmin=xmin, xmax=xmax,
                                      zoom_xmin=7657, zoom_xmax=7663, molecule='O2', plot=True)

        for sp in observations:
            # plt.plot(sp.wave, sp.flux)
            plt.plot(sp.grid, sp.interp_flux)
        plt.show()

    elif example == 2:

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # IMPORT DATA
        file1 = "/HD170740/RED_860/HD170740_w860_redu_20140915_O8.fits"
        file2 = "/HD170740/RED_860/HD170740_w860_redu_20140916_O8.fits"
        file3 = "/HD170740/RED_860/HD170740_w860_redu_20150626_O8.fits"
        file4 = "/HD170740/RED_860/HD170740_w860_redu_20160613_O8.fits"
        file5 = "/HD170740/RED_860/HD170740_w860_redu_20170705_O8.fits"

        xmin = 9575
        xmax = 9580

        sp1 = EdiblesSpectrum(file1)
        sp2 = EdiblesSpectrum(file2)
        sp3 = EdiblesSpectrum(file3)
        sp4 = EdiblesSpectrum(file4)
        sp5 = EdiblesSpectrum(file5)
        observations = [sp1, sp2, sp3, sp4, sp5]

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # CORRECT TELLURIC SHIFT

        observations = telluric_shift(observations, xmin=xmin, xmax=xmax,
                                      zoom_xmin=9576, zoom_xmax=9579, molecule='H2O')

        for sp in observations:
            # plt.plot(sp.wave, sp.flux)
            plt.plot(sp.grid, sp.interp_flux)
        plt.show()

