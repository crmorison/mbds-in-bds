# 'MBD for a B-D process'

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
    divs = np.zeros(pop_i, dtype = int)
    mbs = np.zeros(pop_i, dtype = int)

    while pop < pop_f:
        cell = np.random.randint(0, pop)
        divs = np.append(divs, divs[cell])
        mbs = np.append(mbs, mbs[cell])
        daughters = [cell, -1]
        kid_rand = np.random.poisson(mut_rate, 2)
        pop += 1
        for i in range(len(daughters)):
            divs[daughters[i]] += 1
            if kid_rand[i] > 0:
                mbs[daughters[i]] += int(kid_rand[i])
    return divs, mbs

def simulate(pop_i, pop_f, mut_rate, rls):
    mbss = np.zeros(rls, dtype = list)
    divss = np.zeros(rls, dtype = list)
    print(rls, 'rls from', pop_i, 'to', pop_f, 'cells')
    for i in range(rls):
        divss[i], mbss[i] = ssa_exp(pop_i, pop_f, mut_rate)
    return divss, mbss

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

    if div == 1:
        plt.xlabel('\# of divisions d')
        plt.ylabel('average \# of cells with d divisions')
    else:
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
    if div == 1:
        tikzplotlib.clean_figure()
        # plt.show()
        tikzplotlib.save('av_dd_plot.tex')
    else:
        pois_mbd_y = poissonify([int(i) for i in dd_for_pois], MUT_RATE)
        pois_mbd_x = range(len(pois_mbd_y))
        a = np.nonzero(pois_mbd_y)[0][0]
        b = np.nonzero(pois_mbd_y)[0][-1]
        plt.plot(pois_mbd_x[a:b], pois_mbd_y[a:b], 'b', label = 'predicted')
        plt.legend(loc = 'upper right')
        # plt.show()
        tikzplotlib.clean_figure()
        tikzplotlib_fix_ncols(fig)
        tikzplotlib.save('av_mbd_plot.tex')
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

###############################################################################################

POP_I = 1
POP_F = 10000
MUT_RATE = 10
RLS = 200

divss, mbss = simulate(POP_I, POP_F, MUT_RATE, RLS)
dd_for_pois = av_mbd_plot(RLS, divss, 1, 0)
av_mbd_plot(RLS, mbss, 0, dd_for_pois)
