import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from edibles.sightline import Sightline
from edibles.utils.edibles_spectrum import EdiblesSpectrum

from sourceseparation.wavelength_corr import correct
from sourceseparation.interpolate import interpolate




class Target():


    def __init__(self, sightlines):

        self.dates = [sightline.date for sightline in sightlines]
        self.v_barys = [sightline.v_bary for sightline in sightlines]




