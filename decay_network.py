import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re

def gen_nuclide_matrix(path):
    frame = pd.read_csv(path)
    names = list(set(frame['name']))
    regex = re.compile('[^0-9]')
    names.sort(key=lambda x : (int(regex.sub('', x)) if x!='Neutron' else 0))

    N = len(names)
    print(str(N) + " nuclides in dataset...")
    M = np.zeros((N, N))
    ground_codes = np.empty(N, dtype=object)
    is_ground_state = np.empty(N, dtype=bool)
    
    # very short-lived states list half life in units of energy (width), so we need these later
    hbar = 6.582119569e-16 # in units of eV*s
    ln2 = np.log(2)

    # convert everything to s
    time_units = {
        'ys':1.0e-24,
        'as':1.0e-18,
        'fs':1.0e-15,
        'ps':1.0e-12,
        'ns':1.0e-9,
        'us':1.0e-6,
        'ms':1.0e-3,
        's':1.0,
        'm':60.0,
        'h':60.0*60.0,
        'd':60.0*60.0*24.0,
        'y':60.0*60.0*24.0*365.24,
        'nan':float('nan')
    }

    # figure out which states are ground states (they are printed before others; I rely on this)
    # (excited states are id'd by their energy and partiy)
    def gencode(row):
        return str(row['spinAndParity']) + str(row['bindingEnergy(keV)'])
    for _ , row in frame.iterrows():
        index = names.index(row['name'])
        if ground_codes[index] == None and float(row['levelEnergy(MeV)']) == 0.0:
            ground_codes[index] = gencode(row)
    for _ , row in frame.iterrows():
        index = names.index(row['name'])
        is_ground_state[index] = (gencode(row) == ground_codes[index])

    # dictionary for finding ground state index from (z, n) tuple
    zn_lookup = {}
    for _ , row in frame.iterrows():
        index = names.index(row['name'])
        zn = (int(row['z']), int(row['n']))
        zn_lookup[zn] = index
    a = set()
    # fill the matrix
    for _ , row in frame.iterrows():
        index = names.index(row['name'])

        if is_ground_state[index]:
            half_life = 0.0
            time_unit = str(row['halflifeUnit'])
            if time_unit=='ev':
                half_life = hbar*ln2/(float(row['halflife']))
            elif time_unit=='kev':
                half_life = hbar*ln2/(float(row['halflife'])*1000.0)
            elif time_unit=='mev':
                half_life = hbar*ln2/(float(row['halflife'])*1000000.0)
            elif str(row['halflife'])=='STABLE':
                half_life = np.infty
            else:
                half_life = float(row['halflife'])*time_units[str(row['halflifeUnit'])]

            # diag. elements
            decay_constant = ln2/half_life
            if M[index, index] == 0:
                M[index, index] = -decay_constant

            # off-diag. elements
            # nudat has this stuff in an awkward format: for example, a 100% B- decay probability
            # means a NET 100% beta decay probability, including B-N, B-2N, B-3N, etc. decays
            # so I need to know ALL the probabilities before I can get partial decay rates...
            mode = str(row['decayModes'])

            
            if not ('?' in mode):
                if '=' in mode:
                    a.add(mode)
                    pass
                else:
                    pass
    print(a)
    
    print(ground_codes)
    #plt.imshow(M)
    #plt.show()

# up next: eigenvectors, nice plots

gen_nuclide_matrix("nndc_nudat_data_export_full.csv")
