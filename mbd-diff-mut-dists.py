# 'MBD with different mutational distributions for a B-D process'

# author: Christo Morison

import numpy as np
import matplotlib.pyplot as plt
import decimal as dcm
dcm.getcontext().prec = 500
import math
import tikzplotlib
from pathlib import Path
plt.rc('font', size = 20)

def ssa_exp(pop_i, pop_f, mut_rate):
    pop = pop_i
    mbs = np.zeros(pop_i, dtype = int)
    while pop < pop_f:
        cell = np.random.randint(0, pop)
        mbs = np.append(mbs, mbs[cell])
        daughters = [cell, -1]
        #rando = np.random.randint(0, 2, 2)
        #rands = [mut_rate - 1, mut_rate + 1]
        #kid_rand = [0, 0]
        #for i in range(len(rando)):
        #    kid_rand[i] = rands[rando[i]]
        kid_rand = np.random.geometric(1 / mut_rate, 2)#np.random.poisson(mut_rate, 2)#np.random.randint(1, 2 * mut_rate, 2)#[mut_rate, mut_rate]
        pop += 1
        for i in range(len(daughters)):
            if kid_rand[i] > 0:
                mbs[daughters[i]] += int(kid_rand[i])
    return mbs

def simulate(pop_i, pop_f, mut_rate, rls):
    mbss = np.zeros(rls, dtype = list)
    print(rls, 'rls from', pop_i, 'to', pop_f, 'cells')
    for i in range(rls):
        mbss[i] = ssa_exp(pop_i, pop_f, mut_rate)
    return mbss

def pad_list(list_of_lists, sought_len):
    list_of_lists2 = np.zeros(len(list_of_lists), dtype = list)
    for i in range(len(list_of_lists)):
        list_of_lists2[i] = np.append(list_of_lists[i], np.zeros(sought_len - len(list_of_lists[i])))
    return list_of_lists2

def mbd_prep(mbsy, mbs_hist_maxs):
    max_muts = int(np.max(mbs_hist_maxs))
    burds = np.linspace(0, max_muts, num = max_muts + 1)
    pad_mbs = pad_list(mbsy, max_muts + 1)
    mb_avs = np.average(pad_mbs, axis = 0)
    return burds, mb_avs, pad_mbs

def av_mbd_plot(rls, mbss, div, dd_for_pois):
    mbs_maxs = np.zeros(rls, dtype = int)
    mbsy = np.zeros(rls, dtype = list)
    for i in range(rls):
        mbsy[i] = np.bincount(mbss[i]).astype(np.float64)
        mbs_maxs[i] = mbss[i].max()
    fig = plt.figure(figsize = [8, 6])
    plt.xlabel('\# of mutations j')
    plt.ylabel('average \# of cells with j mutations')
    plt.ticklabel_format(axis = 'both', style = 'sci', scilimits = (-2,3))
    burds, mb_avs, pad_mbs = mbd_prep(mbsy, mbs_maxs)
    sim_mbd_x = [burds[i] for i in np.nonzero(mb_avs)[0]]
    sim_mbd_y = [mb_avs[i] for i in np.nonzero(mb_avs)[0]]
    plt.plot(sim_mbd_x, sim_mbd_y, 'r', label = 'simulated')
    for i in range(0, rls, 7):
        burdsi = [burds[j] for j in np.nonzero(pad_mbs[i])[0]]
        mbsi = [pad_mbs[i][j] for j in np.nonzero(pad_mbs[i])[0]]
        plt.plot(burdsi, mbsi, 'r', alpha = 0.1)
    hist_sim_mbd_x = range(0, int(POP_F - POP_I) * MUT_RATE, MUT_RATE)
    hist_sim_mbd_y = hist_dat(mb_avs, MUT_RATE)
    start1 = np.nonzero(hist_sim_mbd_y)[0][0]
    stop1 =  np.nonzero(hist_sim_mbd_y)[0][-1]
    #plt.plot(hist_sim_mbd_x[start1:stop1], hist_sim_mbd_y[start1:stop1], 'go', label = 'hist')
    pois_mbd_y = poissonify([int(i) for i in hist_sim_mbd_y], MUT_RATE)
    pois_mbd_x = range(len(pois_mbd_y))
    a = np.nonzero(pois_mbd_y)[0][0]
    b = np.nonzero(pois_mbd_y)[0][-1]
    plt.plot(pois_mbd_x[a:b], pois_mbd_y[a:b], 'b', label = 'predicted')
    plt.legend(loc = 'upper right')
    tikzplotlib.clean_figure()
    tikzplotlib_fix_ncols(fig)
    tikzplotlib.save('0av_mbd_plot.tex')
    plt.close()
    return mb_avs

def pad_list(list_of_lists, sought_len):
    list_of_lists2 = np.zeros(len(list_of_lists), dtype = list)
    for i in range(len(list_of_lists)):
        list_of_lists2[i] = np.append(list_of_lists[i], np.zeros(sought_len - len(list_of_lists[i]), dtype = int))
    return list_of_lists2

def pois(k, l):
    return dcm.Decimal(l ** k) * dcm.Decimal(np.exp(-l)) / dcm.Decimal(math.factorial(k))

def poissonify(y_data, mean):
    y_data_pois = np.zeros(len(y_data), dtype = list)
    for i in range(len(y_data)):
        y_data_pois[i] = np.zeros(int((len(y_data) + 2) * mean))
        for k in range(len(y_data_pois[i])):
            y_data_pois[i][k] = y_data[i] * pois(k, i * mean)
    p_data = np.sum(y_data_pois)
    return p_data

def tikzplotlib_fix_ncols(obj):
# workaround for matplotlib 3.6 renamed legend's _ncol to _ncols, which breaks tikzplotlib
    if hasattr(obj, "_ncols"):
        obj._ncol = obj._ncols
    for child in obj.get_children():
        tikzplotlib_fix_ncols(child)

def hist_dat(data_to_bin, bin_size):
    hist = np.zeros(int(len(data_to_bin) / bin_size) + 10)
    for j in range(len(data_to_bin)):
        hist[int(np.rint(np.nextafter(j / bin_size, j / bin_size + 1)))] += data_to_bin[j]
    return hist

###############################################################################################

POP_I = 1
POP_F = 10000
MUT_RATE = 10
RLS = 200

#Path('/Users/christo/Desktop/simplest-mbd-v2/fin=' + str(POP_F) + '/' +
#                str(RLS) + 'rls-' + str(POP_I) + 'ini-' + str(MUT_RATE) + 'mu').mkdir(parents = True, exist_ok = True)


mbss = simulate(POP_I, POP_F, MUT_RATE, RLS)
av_mbd_plot(RLS, mbss, 0, 0)
