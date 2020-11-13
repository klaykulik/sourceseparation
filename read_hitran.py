import astropy.units as u
from specutils.utils.wcs_utils import vac_to_air


def read_hitran(file, molecule='O2'):
    '''This function will read the contents of a HITRAN data file,
    and output wanted parameters.

    Args:
        file (str): The name of the file to be read

    Returns:
        list: A list of dicts, one dict for each line

    Note:
        Currently, this function only works for O2 lines

    '''

    if molecule == 'O2':
        prefix = ' 7'
    if molecule == 'H2O':
        prefix = ' 1'



    pars_list = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith(prefix):

                pars_list.append({'WN': float(line[3:15]),
                                  'S': float(line[16:25]),
                                  'A': float(line[26:35])
                                  })

    return pars_list


def convert(pars_list, air=True):

    new_pars_list = []

    for line_pars in pars_list:
        lam_0 = 1e8 * 1 / line_pars['WN']
        gamma = 1 / line_pars['A']


        # d = gamma * lam_0 * 1e-10 / cst.c.value
        # print(d)



        if line_pars['S'] > 1e-26:
            tau_0 = 0.6
        elif line_pars['S'] > 2e-28:
            tau_0 = 0.03
        elif line_pars['S'] > 1e-29:
            tau_0 = 0.005


        if air:
            air_wave = vac_to_air(lam_0 * u.AA, method='Ciddor1996')
            lam_0 = air_wave.value


        new_pars_list.append({'lam_0': lam_0,
                              'd': gamma,
                              'tau_0': tau_0
                              })


        # print({'lam_0': lam_0, 'd': gamma, 'tau_0': tau_0})
        # print(line_pars['S'])

    return new_pars_list


if __name__ == '__main__':
    pars_list = read_hitran('telluric_lines_HITRAN.txt', molecule='O2')

    # print(pars_list)

    pars_list = convert(pars_list)




    xmin = 7661.75
    xmax = 7669
    # Set tau_cutoff (0.4, 0.02, or 0.0)
    tau_cutoff = 0.02

    # Create linelist
    linelist = []
    for pars in pars_list:
        if (xmin < pars["lam_0"]) and (pars["lam_0"] < xmax):
            if pars["tau_0"] > tau_cutoff:
                linelist.append(pars)
    linelist.reverse()
    print(linelist)

    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"


    from edibles.utils.edibles_spectrum import EdiblesSpectrum
    import matplotlib.pyplot as plt


    sp1 = EdiblesSpectrum(FILE1)
    sp1.getSpectrum(xmin=xmin, xmax=xmax)


    plt.plot(sp1.wave, sp1.flux)


    plt.vlines([line['lam_0'] for line in linelist], 0, 1000)



    plt.show()
