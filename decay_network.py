import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import re
from scipy import sparse
import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
import scipy.linalg
import scipy.sparse.linalg
from scipy.sparse.linalg import eigs
from numpy.linalg import eig

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
    mode_out = np.zeros_like(M)
    legend = {'SF':len(data['nucs']), 'G':(len(data['nucs'])+1)}
    _j = 0
    for x in data['nucs']:
        N = int(x['n'])
        Z = int(x['z'])
        legend[str(Z) + ' ' + str(N)] = _j
        _j += 1
    decay_constant = 0.0
    half_lives = []
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
                    if b > 0.01 and not (this_ID == decayed_ID):
                        M[this_ID][decayed_ID] += decay_constant * b/100.0
                        mode_out[this_ID][decayed_ID] = i
        else:
            decay_constant = np.inf
        half_lives.append(ln2/decay_constant)
    names = []
    for x in data['nucs']:
        N = int(x['n'])
        Z = int(x['z'])
        A = Z+N
        names.append(str(A) + data['elements'][Z])

    return M, names, half_lives, mode_out

def save_matrix(M, path):
    pd.DataFrame(M).to_csv(path, index=False, index_label = False, header=False)

def plot_matrix(M):
    M = np.log10(M)
    M = np.nan_to_num(M, nan=-100.0, posinf=-100.0, neginf=-100.0)
    plt.imshow(M, vmin=-10.0, vmax=10.0, interpolation='none', cmap='magma')
    plt.colorbar()
    plt.show()
    plt.imsave('log10gamma.png',M,cmap='magma', vmin=-10.0, vmax=10.0)

def get_eigen(M):
    w, V = np.linalg.eig(M)
    save_matrix(w, 'evalues.csv')
    save_matrix(V, 'evectors.csv')
    return w, V

def get_coeffs(V, state):
    return np.linalg.solve(V.T, state)

def get_amounts(time, coeffs, eigenvalues, eigenvectors):
    sum = np.zeros_like(eigenvectors[0])
    for (coeff, decay, vec) in zip(coeffs, eigenvalues, eigenvectors):
        print(coeff, time, decay, vec)
        #if coeff != 0.0:
        #    print(time, coeff, decay)
        sum += coeff*np.exp(time*decay)*vec
    return sum

M, names, half_lives, mtb = NUDAT_TO_MATRIX('NUDAT2020.json')
M = np.nan_to_num(M, nan=1e-30, posinf=1e30, neginf=1e30)
sM = scipy.sparse.csc_array(M)
#M = np.nan_to_num(M, nan=np.inf)
save_matrix(mtb, 'modes.csv')

exit()
# code below has numerical errors; need to fix

state_vec = np.zeros((M.shape[0], 1))
state_vec[names.index('14Be')] = 1.0

t = np.power(10.0, np.arange(-10.0, 10.0, 0.2))
y = []
iii = 0
for x in t:
    try:
        ty = np.matmul(scipy.sparse.linalg.expm(x*sM).toarray(), np.array(state_vec).T)
        y.append(ty)
        print(ty)
    except OverflowError:
        break
    print(iii)
    iii += 5
    if iii > 19: break

t = t[:len(y)]


indices = set()
for vector in y:
    j = 0
    for amount in vector:
        if amount > 0.01:
            indices.add(j)
        j += 1
indices = list(indices)
y = (np.array(y).T[indices]).T
names = np.array(names)[indices]

print(t, y)
objs = plt.plot(t, y)
plt.legend(objs, names)
plt.xscale('log',base=10)
plt.show()






#save_matrix(M, "M.csv")
#save_matrix(half_lives, "half_lives.csv")
#print(names)
#get_eigen(M)

# little practice demo thiny
MM = np.array([[-4, 0], [4, 0]])
MM = np.array([[-4, 0, 0, 0, 0], [4, -3, 0, 0, 0], [0, 3, -6, 0, 0], [0, 0, 6, -2, 0], [0, 0, 0, 2, 0]])
w, V = get_eigen(MM)
print(V)
C = get_coeffs(V, np.array([1, 0, 0, 0, 0]).T)
print(C)
t = np.arange(0.0, 2.0, 0.01)
y = [np.matmul(scipy.linalg.expm(x*MM), np.array([1, 0, 0, 0, 0]).T) for x in t]
plt.plot(t, np.array(y).T[0], color='#66ff99')
plt.plot(t, np.array(y).T[1], color='#00b0f0')
plt.plot(t, np.array(y).T[2], color='#ffee44')
plt.plot(t, np.array(y).T[3], color='#ff77bb')
plt.plot(t, np.array(y).T[4], color='#aa66ee')


ax = plt.gca()
ax.set_facecolor('black')
plt.show()

exit()