#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


npts = 250

def trace(g, tr_g):
    tr_g.zero()
    indices = [(b, i, j) for b, i , j in g.all_indices if i == j]
    for ind in indices:
        tr_g += g[ind[0]][ind[1], ind[2]] / len(indices)

pp = PdfPages('g_qp.pdf')
plt.clf()
arch = HDFArchive('measure_qp.h5', 'r')

arch_g = arch['G_qp_tau']
for i_qp, qp in arch_g.items():
    g = arch_g[i_qp]
    g = BlockGf(name_block_generator = [(s, rebinning_tau(b, npts)) for s, b in g])
    npts = len([x for x in g.mesh])
    tr_g = GfImTime(indices = range(1), beta = g.beta, n_points = npts)
    trace(g, tr_g)
    oplot(tr_g, name=i_qp, mode = 'R')
    print i_qp, tr_g.tail[1][0, 0]
g = arch['G_tau']
g = BlockGf(name_block_generator = [(s, rebinning_tau(b, npts)) for s, b in g])
trace(g, tr_g)
oplot(tr_g, name='tot', mode = 'R')
print 'tot', tr_g.tail[1][0, 0]
plt.gca().legend(loc = 'lower center', fontsize = 10)
plt.gca().set_ylim(-1,0.05)
pp.savefig(plt.gcf())
pp.close()
