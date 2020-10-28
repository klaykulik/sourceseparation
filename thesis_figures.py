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


if __name__ == "__main__":

    # reference_frames()

    line_examples()
