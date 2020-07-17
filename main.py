import numpy as np
import matplotlib.pyplot as plt
from specutils.utils.wcs_utils import vac_to_air
from astropy import units as u

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline

from sourceseparation.wavelength_corr import correct
from sourceseparation.read_hitran import read_hitran, convert

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


sightline1 = Sightline(sp1, n_anchors=10)
sightline2 = Sightline(sp2)
sightline3 = Sightline(sp3)
sightline4 = Sightline(sp4)
sightline5 = Sightline(sp5)

sightlines = [sightline1, sightline2, sightline3, sightline4, sightline5]


out = sightline1.model.eval(data=sp1.flux, params=sightline1.model_pars, x=sp1.wave)
resid = sp1.flux - out

plt.plot(sp1.wave, sp1.flux)
plt.plot(sp1.wave, out)
plt.plot(sp1.wave, resid)
plt.show()


# O_2 data from HITRAN
pars_list = convert(read_hitran('telluric_lines_HITRAN.txt'))

# for pars in pars_list:
#     print(pars)

print(pars_list[0])

linelist = []
for pars in pars_list:
    if (xmin < pars['lam_0']) and (pars['lam_0'] < xmax):

        linelist.append(pars)






plt.plot(sp1.wave, sp1.flux)


tau_cutoff = 0.00



for pars in linelist:

    if pars['tau_0'] > tau_cutoff:
        print(pars)
        plt.scatter(pars['lam_0'], 100)

plt.show()


sightline1.add_source(name='O2', similar={'b': 0.6})

for i in range(len(linelist)):
    line = linelist[i]
    if line['tau_0'] > tau_cutoff:
        name = 'line' + str(i)
        sightline1.add_line(name=name, source='O2', pars=line)
        # print(sightline1.model_pars)
        par_name_d = 'O2_' + name + '_d'
        sightline1.model_pars[par_name_d].set(value=0.05)

out = sightline1.model.eval(data=sp1.flux, params=sightline1.model_pars, x=sp1.wave)
resid = sp1.flux - out

plt.plot(sp1.wave, sp1.flux)
plt.plot(sp1.wave, out)
plt.plot(sp1.wave, resid)
plt.show()


sightline1.fit(report=True, plot=True)



out = sightline1.model.eval(data=sp1.flux, params=sightline1.result.params, x=sp1.wave)
resid = sp1.flux - out

plt.plot(sp1.wave, sp1.flux)
plt.plot(sp1.wave, out)
plt.plot(sp1.wave, resid)



for pars in linelist:

    if pars['tau_0'] > tau_cutoff:

        if pars['tau_0'] == 0.5:
            plt.scatter(pars['lam_0'], 100, marker='.')

        else:
            plt.scatter(pars['lam_0'], 1500, marker='.')

plt.show()
