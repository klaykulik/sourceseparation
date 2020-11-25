from __future__ import print_function
import numpy as np
from scipy.signal import find_peaks

from edibles.models.create_model import createCont
from edibles.models.model import Sightline
from edibles.fitter import fit


def normalize(data, star_name, year):

    prominence = (np.max(data[1]) - np.min(data[1])) * 0.1
    peaks, _ = find_peaks(-data[1], prominence=prominence)
    peak_wavelengths = data[0].iloc[peaks]

    cont = createCont(data, n_points=4)

    sightline = Sightline(star_name=star_name, cont=cont)
    for i in range(len(peak_wavelengths)):
        source_name = "Source " + str(i + 1)
        line_name = "line" + str(i + 1)

    # sightline.addSource(source_name='source_name1', b=3.0, d=3)
    # sightline.addLine(name='line_name1', lam_0=6285, tau_0=0.3)

    fit_m, params = fit(star_name, data, sightline.model, breakdown=True, silent=True)

    f_norm = data[1] / cont(data[0])
    norm_data = [data[0], f_norm]

    return norm_data
