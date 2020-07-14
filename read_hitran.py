
def read_hitran(file):


    # lines = {}


    with open(file, 'r') as f:
        for line in f:
            if line.startswith(' 71'):
                wave = line[3:16]
                print(wave)
                S = line[16:25]
                print(S)
                A = line[26:35]
                print(A)
                Gamma = 1/float(A)
                print(Gamma)
                # lines.append


read_hitran('telluric_lines_HITRAN.txt')

a = '2.2e-01'
b = float(a)
print(b)