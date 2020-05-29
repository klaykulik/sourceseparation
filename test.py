import numpy as np
import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models.create_model import createCont
from edibles.models.model import Sightline
from edibles.fitter import fit

from wavelength_corr import correctWavelength


FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
FILE2 = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
FILE3 = "/HD170740/RED_860/HD170740_w860_redl_20150626_O12.fits"
FILE4 = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
FILE5 = "/HD170740/RED_860/HD170740_w860_redl_20170705_O12.fits"

sp1 = EdiblesSpectrum(FILE1)
subset1 = sp1.getSpectrum(xmin=7661.5, xmax=7669)
# subset1.flux = subset1.flux / np.max(subset1.flux)
data = (subset1.wave, subset1.flux)
cont = createCont(data, n_points=4)
sightline1 = Sightline(star_name=sp1.target, cont=cont)
sightline1.addSource(source_name="Telluric", b=1.1, d=0.05)
sightline1.addLine(name="tell1", lam_0=7664.8, tau_0=0.6)
sightline1.addLine(name="tell2", lam_0=7666, tau_0=0.6)
sightline1.addSource(source_name="Interstellar", b=1.1, d=0.02)
sightline1.addLine(name="int1", lam_0=7665.2, tau_0=0.2)
sightline1.addLine(name="int2", lam_0=7665.3, tau_0=0.01)
fit_model = fit(sp1.target, data, sightline1.model, breakdown=False, silent=True)
line_vals1 = np.array([line.lam_0.val for line in sightline1.lines["all"]])
y1 = np.ones_like(line_vals1) * sp1.v_bary

sp2 = EdiblesSpectrum(FILE2)
subset2 = sp2.getSpectrum(xmin=7661.5, xmax=7669)
# subset2.flux = subset2.flux / np.max(subset2.flux)
data = (subset2.wave, subset2.flux)
cont = createCont(data, n_points=4)
sightline2 = Sightline(star_name=sp2.target, cont=cont)
sightline2.addSource(source_name="Telluric", b=1.1, d=0.05)
sightline2.addLine(name="tell1", lam_0=7664.8, tau_0=0.6)
sightline2.addLine(name="tell2", lam_0=7666, tau_0=0.6)
sightline2.addSource(source_name="Interstellar", b=1.1, d=0.02)
sightline2.addLine(name="int1", lam_0=7665.3, tau_0=0.2)
sightline2.addLine(name="int2", lam_0=7665.4, tau_0=0.01)
fit_model = fit(sp2.target, data, sightline2.model, breakdown=False, silent=True)
line_vals2 = np.array([line.lam_0.val for line in sightline2.lines["all"]])
y2 = np.ones_like(line_vals2) * sp2.v_bary

sp3 = EdiblesSpectrum(FILE3)
subset3 = sp3.getSpectrum(xmin=7661.5, xmax=7669)
# subset3.flux = subset3.flux / np.max(subset3.flux)
data = (subset3.wave, subset3.flux)
cont = createCont(data, n_points=4)
sightline3 = Sightline(star_name=sp3.target, cont=cont)
sightline3.addSource(source_name="Telluric", b=1.1, d=0.05)
sightline3.addLine(name="tell1", lam_0=7664.8, tau_0=0.6)
sightline3.addLine(name="tell2", lam_0=7666, tau_0=0.6)
sightline3.addSource(source_name="Interstellar", b=1.1, d=0.02)
sightline3.addLine(name="int1", lam_0=7664.5, tau_0=0.2)
sightline3.addLine(name="int2", lam_0=7664.6, tau_0=0.01)
fit_model = fit(sp3.target, data, sightline3.model, breakdown=False, silent=True)
line_vals3 = np.array([line.lam_0.val for line in sightline3.lines["all"]])
y3 = np.ones_like(line_vals3) * sp3.v_bary

sp4 = EdiblesSpectrum(FILE4)
subset4 = sp4.getSpectrum(xmin=7661.5, xmax=7669)
# subset4.flux = subset4.flux / np.max(subset4.flux)
data = (subset4.wave, subset4.flux)
cont = createCont(data, n_points=4)
sightline4 = Sightline(star_name=sp4.target, cont=cont)
sightline4.addSource(source_name="Telluric", b=1.1, d=0.05)
sightline4.addLine(name="tell1", lam_0=7664.8, tau_0=0.6)
sightline4.addLine(name="tell2", lam_0=7666, tau_0=0.6)
sightline4.addSource(source_name="Interstellar", b=1.1, d=0.02)
sightline4.addLine(name="int1", lam_0=7664.4, tau_0=0.2)
sightline4.addLine(name="int2", lam_0=7664.5, tau_0=0.01)
fit_model = fit(sp4.target, data, sightline4.model, breakdown=False, silent=True)
line_vals4 = np.array([line.lam_0.val for line in sightline4.lines["all"]])
y4 = np.ones_like(line_vals4) * sp4.v_bary


sp5 = EdiblesSpectrum(FILE5)
subset5 = sp5.getSpectrum(xmin=7661.5, xmax=7669)
# subset5.flux = subset5.flux / np.max(subset5.flux)
data = (subset5.wave, subset5.flux)
cont = createCont(data, n_points=4)
sightline5 = Sightline(star_name=sp5.target, cont=cont)
sightline5.addSource(source_name="Telluric", b=1.1, d=0.05)
sightline5.addLine(name="tell1", lam_0=7664.8, tau_0=0.6)
sightline5.addLine(name="tell2", lam_0=7666, tau_0=0.6)
sightline5.addSource(source_name="Interstellar", b=1.1, d=0.02)
sightline5.addLine(name="int1", lam_0=7664.6, tau_0=0.2)
sightline5.addLine(name="int2", lam_0=7664.65, tau_0=0.01)
fit_model = fit(sp5.target, data, sightline5.model, breakdown=False, silent=True)
line_vals5 = np.array([line.lam_0.val for line in sightline5.lines["all"]])
y5 = np.ones_like(line_vals5) * sp5.v_bary


plt.scatter(line_vals1, y1)
plt.scatter(line_vals2, y2)
plt.scatter(line_vals3, y3)
plt.scatter(line_vals4, y4)
plt.scatter(line_vals5, y5)
plt.show()

# print(i for i in line_vals1)
# x = []
# x.append(i for i in line_vals1)
# x.append([i for i in line_vals2])
# x.append([i for i in line_vals3])
# x.append([i for i in line_vals4])
# x.append([i for i in line_vals5])

# x = np.array(x)

# print(x.T)


# # print(y1, y2, y3, y4, y5)

# y = []
# y.append([i for i in y1])
# y.append([i for i in y2])
# y.append([i for i in y3])
# y.append([i for i in y4])
# y.append([i for i in y5])
# y = np.array(y)
# print(y.T)
