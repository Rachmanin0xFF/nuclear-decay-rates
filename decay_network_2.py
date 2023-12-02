import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import re

import networkx as nx

def show_graph_with_labels(adjacency_matrix):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=500, with_labels=False)
    plt.show()

def NUDAT_TO_MATRIX(path):
    units = set()
    ln2 = np.log(2)
    time_units = {
        'ys':1.0e-24,
        'zs':1.0e-21,
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
        'ky':60.0*60.0*24.0*365.24*1.0e3,
        'My':60.0*60.0*24.0*365.24*1.0e6,
        'Gy':60.0*60.0*24.0*365.24*1.0e9,
        'Ty':60.0*60.0*24.0*365.24*1.0e12,
        'Py':60.0*60.0*24.0*365.24*1.0e15,
        'Ey':60.0*60.0*24.0*365.24*1.0e18,
        'Zy':60.0*60.0*24.0*365.24*1.0e21,
        'Yy':60.0*60.0*24.0*365.24*1.0e24,
        'nan':float('nan')
    }
    mode_table = pd.read_csv('mode_table.csv')
    data = None
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    if data == None:
        return None

    # last two are SF, G
    M = np.zeros(shape=(len(data['nucs'])+2, len(data['nucs'])+2))
    legend = {'SF':len(data['nucs']), 'G':(len(data['nucs'])+1)}
    _j = 0
    for x in data['nucs']:
        N = int(x['n'])
        Z = int(x['z'])
        legend[str(Z) + ' ' + str(N)] = _j
        _j += 1
    print(legend)
    
    for x in data['nucs']:
        N = int(x['n'])
        Z = int(x['z'])
        A = N+Z
        this_ID = legend[str(Z) + ' ' + str(N)]
        if 'dm' in x.keys():
            #print("----- " + x['z'] + ' ' + x['n'] + ' ------')
            # Branching ratios and uncertainties
            BR = [0] * len(x['dm'])
            BR_error = [0] * len(x['dm'])
            modes = [0] * len(x['dm'])

            ONM_index = -1
            for i in range(len(x['dm'])):
                if 'a' in x['dm'][i].keys():
                    modes[i] = int(x['dm'][i]['a'])
                if 'b' in x['dm'][i].keys():
                    text = x['dm'][i]['b']
                    if 'u003d?' in text:
                        ONM_index = i
                    text = text.replace('u003d?','0.0')
                    text = text.replace('?', '0')
                    text = text.replace('~', '')
                    text = text.replace('LT', '') # not really the best way to do things, but I'm not sure what is
                    text = text.replace('GT', '')
                    text = text.replace('u003d', '')
                    text = text.replace('#', '')
                    ratio = 0.0
                    error = 0.0
                    try:
                        ratio = float(text.split(' ')[0])
                    except ValueError:
                        pass
                    try:
                        error = float(text.split(' ')[2])
                    except:
                        pass

                    BR[i] = ratio
                    BR_error[i] = error
            BR = np.array(BR)
            BR_sum = np.sum(BR)
            if ONM_index != -1:
                BR[ONM_index] = max(0.0, 100.0 - BR_sum)
            BR_sum = np.sum(BR)
            if BR_sum != 100:
                if BR_sum > 100.0 and all(b < 100.0 for b in BR):
                    BR /= 0.01*BR_sum
                else:
                    BR_sum = np.sum(BR)
                    if BR_sum == 0.0:
                        BR[0] = 100.0
                    elif BR_sum < 100.0 and BR[0] == 0:
                        BR[0] = 100.0 - BR_sum
                    elif BR_sum < 100.0 and BR[0] != 0 and len(BR) > 1 and BR[1] == 0:
                        BR[1] = 100.0 - BR_sum
                    elif BR_sum > 100.0:
                        BR[0] = 200.0 - BR_sum
                
                BR_sum =  np.sum(BR)
            
            h_text = x['h']
            h_text = h_text.replace('#', '')
            h_text = h_text.replace('~', '')
            h_text = h_text.replace('lt', '')
            h_text = h_text.replace('gt', '')

            decay_constant = 0.0
            
            if 'stable' not in h_text:
                value = h_text.split(' ')[0]
                #print(h_text)
                unit = h_text.split(' ')[1]
                units.add(unit)
                if value=='' or value=='p-unst':
                    value = 0.0
                else:
                    value = float(value)
                if unit == '':
                    half_life = 0.0
                else:
                    half_life = float(value)*time_units[unit]
                decay_constant = ln2/half_life
            
            M[this_ID][this_ID] = -decay_constant
            for i, b in zip(modes, BR):
                mode_name = mode_table['Decay Type'][i]
                dN = mode_table['dN'][i]
                dZ = mode_table['dZ'][i]
                decayed_ID = 0
                if dN == 'SF':
                    decayed_ID = legend['SF']
                else:
                    new_N = N + int(dN)
                    new_Z = Z + int(dZ)
                    try:
                        decayed_ID = legend[str(new_Z) + ' ' + str(new_N)]
                    except KeyError:
                        b = 0 # we don't accept decays to things that don't exist...
                        decayed_ID = 0
                M[this_ID][decayed_ID] = decay_constant * b/100.0
    M[np.isnan(M)] = 0.0
    M[M>0.0] = 1.0
    #G=nx.from_numpy_array(M)
    #nx.draw(G)
            
    
    
    #print(data['nucs'][1050])
    #print(data['decays'])

M = NUDAT_TO_MATRIX('NUDAT2020.json')

