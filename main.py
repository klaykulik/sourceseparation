import numpy as np
import matplotlib.pyplot as plt
import copy

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline

from sourceseparation.read_hitran import read_hitran, convert
from sourceseparation.interpolate import interpolate


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# IMPORT DATA

file1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
file2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
file3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
file4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
file5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

xmin = 7642
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
# plt.show()

# for sightline in sightlines:
    plt.plot(sightline.grid, sightline.interp_bary_flux)
plt.show()


# O_2 data from HITRAN
pars_list = convert(read_hitran('telluric_lines_HITRAN.txt'))

linelist = []
for pars in pars_list:
    if (xmin < pars['lam_0']) and (pars['lam_0'] < xmax):

        linelist.append(pars)


# Set tau_cutoff (0.4, 0.02, or 0.0)
tau_cutoff = 0.4

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TELLURIC

resids = []
outs = []
for sightline in sightlines:

    sightline.add_source(name='O2', similar={'b': 0.6})

    for i in range(len(linelist)):
        line = linelist[i]
        if line['tau_0'] > tau_cutoff:
            name = 'line' + str(i)
            sightline.add_line(name=name, source='O2', pars=line)
            # print(sightline.model_pars)
            par_name_d = 'O2_' + name + '_d'
            sightline.model_pars[par_name_d].set(value=0.05)

    out = sightline.model.eval(data=sightline.interp_flux,
                               params=sightline.model_pars,
                               x=sightline.grid)
    resid = sightline.interp_flux - out

    # PLOT unfit lines
    # plt.plot(sightline.grid, sightline.interp_flux)
    # plt.plot(sightline.grid, out)
    # plt.plot(sightline.grid, resid)
    # plt.show()

    sightline.fit(report=True, plot=False)

    out = sightline.model.eval(data=sightline.interp_flux,
                               params=sightline.result.params,
                               x=sightline.grid)
    resid = sightline.interp_flux - out

    # PLOT fit lines in one sightline
    # plt.plot(sightline.grid, sightline.interp_flux)
    # plt.plot(sightline.grid, out)
    # plt.plot(sightline.grid, resid)



    resids.append(resid)
    outs.append(out)

# PLOT fit lines in all sightlines
for i in range(len(resids)):
    # plt.plot(sightlines[i].grid, outs[i], marker='.')
    plt.plot(sightlines[i].grid, resids[i], marker='.')
    # print(len(sightlines[i].grid), len(resids[i]))
plt.show()


coadd = np.ones_like(resids[0])

for resid in resids:
    coadd *= resid

plt.plot(sightlines[0].grid, coadd)
plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# WHAT NEXT
