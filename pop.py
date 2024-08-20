# 'Expected population for a B-D process'

# author: Christo Morison

import numpy as np
import matplotlib.pyplot as plt
import decimal as dcm
dcm.getcontext().prec = 100
import math
import tikzplotlib
from pathlib import Path
plt.rc('font', size = 20)

def tikzplotlib_fix_ncols(obj):
# workaround for matplotlib 3.6 renamed legend's _ncol to _ncols, which breaks tikzplotlib
    if hasattr(obj, "_ncols"):
        obj._ncol = obj._ncols
    for child in obj.get_children():
        tikzplotlib_fix_ncols(child)

p = np.linspace(0.5, 0.95, num = 30)
N = [100, 200, 400]

extinctprob = np.zeros(len(N), dtype = list)
xi = np.zeros(len(N), dtype = list)
predti1N = np.zeros(len(N), dtype = list)

for NN in range(len(N)):
    predti1N[NN] = np.zeros(len(p), dtype = list)
    extinctprob[NN] = np.zeros((len(p), N[NN] + 1))
    xi[NN] = np.zeros((len(p), len(extinctprob[NN][0])))

    for pp in range(len(p)):
        for i in range(1, len(extinctprob[NN][pp])):
            extinctprob[NN][pp][i] = np.sum([(dcm.Decimal(math.factorial(2 * n)) / (dcm.Decimal(math.factorial(n)) *dcm.Decimal( math.factorial(n + 1)))) * dcm.Decimal(p[pp] ** n) * dcm.Decimal((1 - p[pp]) ** (n + 1)) for n in range(i)])
        for i in range(len(xi[NN][pp])):
            xi[NN][pp][i] = 1 - extinctprob[NN][pp][i]

        predti1N[NN][pp] = [extinctprob[NN][pp][N[NN]]]
        for n in range(1, N[NN], 2):
            predti1N[NN][pp].append(dcm.Decimal(n / (N[NN] + 1)) * (dcm.Decimal(math.factorial(N[NN] + 1)) / (dcm.Decimal(math.factorial((N[NN] - (n - 1)) / 2)) * dcm.Decimal(math.factorial((N[NN] + (n + 1)) / 2))) * dcm.Decimal(p[pp] ** ((N[NN] + (n - 1)) / 2)) * dcm.Decimal((1 - p[pp]) ** ((N[NN] - (n - 1)) / 2))))
            predti1N[NN][pp].append(0)

fig = plt.figure(figsize = [8, 6])
for NN in range(len(N)):
    plt.ylabel('E[N_i | N_i > 0]')
    plt.xlabel('b')
    plt.plot(p, [np.sum([dcm.Decimal(n) * dcm.Decimal(predti1N[NN][pp][n]) for n in range(len(predti1N[NN][pp]))]) / dcm.Decimal(xi[NN][pp][N[NN]]) for pp in range(len(p))], ':')
    plt.plot(p, [1 + (2 * p[pp] - 1) * N[NN] for pp in range(len(p))], label = str(N[NN]))
plt.legend(loc = 'upper left')

tikzplotlib.clean_figure()
tikzplotlib_fix_ncols(fig)
tikzplotlib.save('eni.tex')
plt.close()
