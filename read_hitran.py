import astropy.units as u
from specutils.utils.wcs_utils import vac_to_air


def read_hitran(file):
    '''This function will read the contents of a HITRAN data file,
    and output wanted parameters.

    Args:
        file (str): The name of the file to be read

    Returns:
        list: A list of dicts, one dict for each line

    Note:
        Currently, this function only works for O2 lines

    '''


    pars_list = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith(' 7'):

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
            tau_0 = 0.5
        elif line_pars['S'] > 1e-27:
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




    return new_pars_list


if __name__ == '__main__':
    pars_list = read_hitran('telluric_lines_HITRAN.txt')

    print(pars_list[0])

    new_list = convert(pars_list)

    for pars in new_list:
        print(pars)
