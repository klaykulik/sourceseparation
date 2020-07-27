import numpy as np
import matplotlib.pyplot as plt
import copy
import astropy.constants as cst

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


tell_sightlines = copy.deepcopy(sightlines)
bary_sightlines = copy.deepcopy(sightlines)

# for sightline in sightlines:
#     plt.plot(sightline.grid, sightline.interp_flux)
# # plt.show()

# # for sightline in sightlines:
#     plt.plot(sightline.grid, sightline.interp_bary_flux)
# plt.show()


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
print(linelist)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TELLURIC

resids = []
outs = []
tell_errs = []
tell_lams = []
for sightline in tell_sightlines:

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

    sightline.fit(report=True, plot=True)

    out = sightline.model.eval(data=sightline.interp_flux,
                               params=sightline.result.params,
                               x=sightline.grid)
    resid = sightline.interp_flux - out


    # PLOT fit lines in one sightline
    # plt.plot(sightline.grid, sightline.interp_flux)
    # plt.plot(sightline.grid, out)
    # plt.plot(sightline.grid, resid)



    errs = []
    lams = []
    for name in sightline.model.param_names:
        if name[-5:] == 'lam_0':
            errs.append(sightline.result.params[name].stderr)
            lams.append(sightline.result.params[name].value)

    tell_errs.append(errs)
    tell_lams.append(lams)



    resids.append(resid)
    outs.append(out)


# # PLOT fit lines in all sightlines
for i in range(len(resids)):
    # plt.plot(tell_sightlines[i].grid, tell_sightlines[i].interp_flux, marker='.')
    # plt.plot(tell_sightlines[i].grid, outs[i], marker='.')
    plt.plot(tell_sightlines[i].grid, resids[i], marker='.')

plt.show()


coadd = np.ones_like(resids[0])

for resid in resids:
    coadd *= resid

plt.plot(tell_sightlines[0].grid, coadd)
plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # WHAT NEXT

# for each sightline:
#     tell_data, bary_data, tell_model, tell_out, tell_resid

# we also have:
#     coadded


# try to get bary_resid:
#     bary_resid = bary_data - bary_out
#     bary_out = fit(model with bary_data)



# tell_sightlines[0].add_source(name='O2', similar={'b': 0.6})

# for i in range(len(linelist)):
#     line = linelist[i]
#     if line['tau_0'] > tau_cutoff:
#         name = 'line' + str(i)
#         tell_sightlines[0].add_line(name=name, source='O2', pars=line)
#         # print(tell_sightlines[0].model_pars)
#         par_name_d = 'O2_' + name + '_d'
#         tell_sightlines[0].model_pars[par_name_d].set(value=0.05)

# out = tell_sightlines[0].model.eval(data=tell_sightlines[0].interp_flux,
#                                     params=tell_sightlines[0].model_pars,
#                                     x=tell_sightlines[0].grid)
# resid = tell_sightlines[0].interp_flux - out

# # PLOT unfit lines
# # plt.plot(tell_sightlines[0].grid, tell_sightlines[0].interp_flux)
# # plt.plot(tell_sightlines[0].grid, out)
# # plt.plot(tell_sightlines[0].grid, resid)
# # plt.show()

# result = tell_sightlines[0].fit(report=True, plot=False)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# plt.figure()
# for i in range(len(linelist)):
#     line = linelist[i]
#     if line['tau_0'] > tau_cutoff:
#         plt.scatter(line['lam_0'], 0, c='k')
#         plt.vlines(x=line['lam_0'], ymin=-30, ymax=10)

# # print(sightline)


# for i in range(len(tell_sightlines)):
#     sightline = tell_sightlines[i]

#     for name in sightline.model.param_names:
#         if name[-5:] == 'lam_0':

#             plt.scatter(sightline.result.params[name].value, sightline.v_bary)


# plt.show()



reals = [linelist[i]['lam_0'] for i in range(len(linelist))]



print(reals)



cs = ['C0', 'C1', 'C2', 'C3', 'C4']
for i in range(len(tell_sightlines)):
    lams = tell_lams[i]
    errs = tell_errs[i]

    shift_resids = np.subtract(reals, lams)


    plt.errorbar(reals, shift_resids, yerr=errs, fmt='o', c=cs[i], label=tell_sightlines[i].datetime.date())

plt.scatter(reals, np.zeros_like(reals), c='k', label='HITRAN data')

plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Residual fit difference (Angstroms)')
plt.title("HD170740")
plt.legend()
plt.show()





