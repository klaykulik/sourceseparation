import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


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
prominence1 = (np.max(subset1.flux) - np.min(subset1.flux)) * 0.1
peaks1, _ = find_peaks(-subset1.flux, prominence=prominence1)
peak_wavelengths1 = [subset1.wave.iloc[i] for i in peaks1]
peak_fluxes1 = [subset1.flux.iloc[i] for i in peaks1]


sp2 = EdiblesSpectrum(FILE2)
subset2 = sp2.getSpectrum(xmin=7661.5, xmax=7669)
# subset2.flux = subset2.flux / np.max(subset2.flux)
prominence2 = (np.max(subset2.flux) - np.min(subset2.flux)) * 0.1
peaks2, _ = find_peaks(-subset2.flux, prominence=prominence2)
peak_wavelengths2 = [subset2.wave.iloc[i] for i in peaks2]
peak_fluxes2 = [subset2.flux.iloc[i] for i in peaks2]


sp3 = EdiblesSpectrum(FILE3)
subset3 = sp3.getSpectrum(xmin=7661.5, xmax=7669)
# subset3.flux = subset3.flux / np.max(subset3.flux)
prominence3 = (np.max(subset3.flux) - np.min(subset3.flux)) * 0.1
peaks3, _ = find_peaks(-subset3.flux, prominence=prominence3)
peak_wavelengths3 = [subset3.wave.iloc[i] for i in peaks3]
peak_fluxes3 = [subset3.flux.iloc[i] for i in peaks3]


sp4 = EdiblesSpectrum(FILE4)
subset4 = sp4.getSpectrum(xmin=7661.5, xmax=7669)
# subset4.flux = subset4.flux / np.max(subset4.flux)
prominence4 = (np.max(subset4.flux) - np.min(subset4.flux)) * 0.1
peaks4, _ = find_peaks(-subset4.flux, prominence=prominence4)
peak_wavelengths4 = [subset4.wave.iloc[i] for i in peaks4]
peak_fluxes4 = [subset4.flux.iloc[i] for i in peaks4]


sp5 = EdiblesSpectrum(FILE5)
subset5 = sp5.getSpectrum(xmin=7661.5, xmax=7669)
# subset5.flux = subset5.flux / np.max(subset5.flux)
prominence5 = (np.max(subset5.flux) - np.min(subset5.flux)) * 0.1
peaks5, _ = find_peaks(-subset5.flux, prominence=prominence5)
peak_wavelengths5 = [subset5.wave.iloc[i] for i in peaks5]
peak_fluxes5 = [subset5.flux.iloc[i] for i in peaks5]


plt.plot(subset1.wave, subset1.flux)
plt.scatter(peak_wavelengths1, peak_fluxes1)

plt.plot(subset2.wave, subset2.flux)
plt.scatter(peak_wavelengths2, peak_fluxes2)

plt.plot(subset3.wave, subset3.flux)
plt.scatter(peak_wavelengths3, peak_fluxes3)

plt.plot(subset4.wave, subset4.flux)
plt.scatter(peak_wavelengths4, peak_fluxes4)

plt.plot(subset5.wave, subset5.flux)
plt.scatter(peak_wavelengths5, peak_fluxes5)

plt.show()


y_vlist1 = np.ones_like(peak_wavelengths1) * sp1.v_bary
y_vlist2 = np.ones_like(peak_wavelengths2) * sp2.v_bary
y_vlist3 = np.ones_like(peak_wavelengths3) * sp3.v_bary
y_vlist4 = np.ones_like(peak_wavelengths4) * sp4.v_bary
y_vlist5 = np.ones_like(peak_wavelengths5) * sp5.v_bary


plt.scatter(peak_wavelengths1, y_vlist1, label=sp1.date[0:10])
plt.scatter(peak_wavelengths2, y_vlist2, label=sp2.date[0:10])
plt.scatter(peak_wavelengths3, y_vlist3, label=sp3.date[0:10])
plt.scatter(peak_wavelengths4, y_vlist4, label=sp4.date[0:10])
plt.scatter(peak_wavelengths5, y_vlist5, label=sp5.date[0:10])

plt.xlabel('Peak Wavelength')
plt.ylabel('Barycentric Velocity')
plt.legend()


plt.axvline(x=peak_wavelengths2[0])
plt.axvline(x=peak_wavelengths2[-1])


linex = [peak_wavelengths1[1], peak_wavelengths2[1], peak_wavelengths3[0], peak_wavelengths4[0], peak_wavelengths5[0]]
liney = [sp1.v_bary, sp2.v_bary, sp3.v_bary, sp4.v_bary, sp5.v_bary]

z = np.polyfit(linex, liney, 1)
p = np.poly1d(z)
plt.plot(linex, p(linex))


plt.show()


datax = peak_wavelengths1 + peak_wavelengths2 + peak_wavelengths3 + peak_wavelengths4 + peak_wavelengths5
datay = list(y_vlist1) + list(y_vlist2) + list(y_vlist3) + list(y_vlist4) + list(y_vlist5)
classification = ['blue', 'red', 'blue', 'blue', 'red', 'blue', 'red', 'blue', 'blue', 'red', 'blue', 'blue', 'red', 'blue', 'blue']
c_num = [0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]

print(datax)
print(datay)

plt.scatter(datax, datay, color=classification)
plt.show()




