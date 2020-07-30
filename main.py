import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline

from sourceseparation.read_hitran import read_hitran, convert


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
sp1.getSpectrum(xmin, xmax)
sp2 = EdiblesSpectrum(file2)
sp2.getSpectrum(xmin, xmax)
sp3 = EdiblesSpectrum(file3)
sp3.getSpectrum(xmin, xmax)
sp4 = EdiblesSpectrum(file4)
sp4.getSpectrum(xmin, xmax)
sp5 = EdiblesSpectrum(file5)
sp5.getSpectrum(xmin, xmax)
observations = [sp1, sp2, sp3, sp4, sp5]

sightline1 = Sightline(sp1, n_anchors=5)
sightline2 = Sightline(sp2, n_anchors=5)
sightline3 = Sightline(sp3, n_anchors=5)
sightline4 = Sightline(sp4, n_anchors=5)
sightline5 = Sightline(sp5, n_anchors=5)
sightlines = [sightline1, sightline2, sightline3, sightline4, sightline5]

for sightline in sightlines:
    plt.plot(sightline.grid, sightline.interp_flux)
    plt.plot(sightline.grid, sightline.interp_bary_flux)
plt.show()

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
# print(linelist)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIT TELLURIC HITRAN LINES
resids = []
outs = []
errs_list = []
lams_list = []
for sightline in sightlines:
    sightline.add_source(name='O2', similar={'b': 0.6})

    for i in range(len(linelist)):
        line = linelist[i]
        if line['tau_0'] > tau_cutoff:
            name = 'line' + str(i)
            sightline.add_line(name=name, source='O2', pars=line)
            par_name_d = 'O2_' + name + '_d'
            sightline.model_pars[par_name_d].set(value=0.05)

    sightline.fit(report=True, plot=False)

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

# PLOT fit lines in all sightlines
for i in range(len(resids)):
    plt.plot(sightlines[i].wave, resids[i], marker='.')

plt.show()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CORRECT WAVELENGTH SHIFT FROM FIT
reals = [linelist[i]['lam_0'] for i in range(len(linelist))]
print(reals)
shift_resids = [np.subtract(reals, lams) for lams in lams_list]

cs = ['C0', 'C1', 'C2', 'C3', 'C4']
for i in range(len(sightlines)):
    lams = lams_list[i]
    errs = errs_list[i]
    shift = shift_resids[i]

    plt.errorbar(reals, shift, yerr=errs, fmt='o', c=cs[i], label=sightlines[i].datetime.date())

plt.scatter(reals, np.zeros_like(reals), c='k', label='HITRAN data')

plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Residual fit difference (Angstroms)')
plt.title("HD170740")
plt.legend()
plt.show()

for i in range(len(sightlines)):

    weights = np.mean(shift_resids[i]) - np.abs(np.mean(shift_resids[i]) - shift_resids[i])
    weights = 1 / shift_resids[i]

    spl = UnivariateSpline(reals, shift_resids[i])
    # spl.set_smoothing_factor(5)

    plt.scatter(reals, shift_resids[i], c=cs[i])
    plt.plot(reals, spl(reals), linestyle='dashed')
    plt.plot(sightlines[i].wave, spl(sightlines[i].wave), c=cs[i])
plt.show()



shifted_grids = []
for i in range(len(sightlines)):

    weights = np.mean(shift_resids[i]) - np.abs(np.mean(shift_resids[i]) - shift_resids[i])
    weights = 1 / shift_resids[i]

    spl = UnivariateSpline(reals, shift_resids[i])
    # spl.set_smoothing_factor(5)

    # plt.plot(sightlines[i].wave + spl(sightlines[i].wave), sightlines[i].flux, c=cs[i])



    shifted_grids.append(sightlines[i].wave + spl(sightlines[i].wave))



for i in range(len(sightlines)):


    plt.plot(shifted_grids[i], sightlines[i].flux, label=sightlines[i].datetime.date())


plt.scatter(reals, np.zeros_like(reals), c='k', label='HITRAN data')
plt.legend()
plt.show()
