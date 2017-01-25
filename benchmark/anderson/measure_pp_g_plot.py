#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np, matplotlib
from mpl_to_latex.matplotlib_to_latex import set_prl_parameters
set_prl_parameters(72, 390)

npts = 500

def trace(g, tr_g):
    tr_g.zero()
    indices = [(b, i, j) for b, i , j in g.all_indices if i == j]
    for ind in indices:
        tr_g += g[ind[0]][ind[1], ind[2]] / len(indices)

pp = PdfPages('g_pp.pdf')
plt.clf()
arch = HDFArchive('measure_pp.h5', 'r')
fullgbygpp = None
arch_g = arch['G_pp_tau']
n_pp = len([k for k in arch_g.keys()])
colors = [matplotlib.cm.jet(i/float(n_pp)) for i in range(n_pp)]
for i_pp, pp in arch_g.items():
    g = arch_g[i_pp]
    if fullgbygpp is None:
        fullgbygpp = g.copy()
    else:
        fullgbygpp += g
    g = BlockGf(name_block_generator = [(s, rebinning_tau(b, npts)) for s, b in g])
    tr_g = GfImTime(indices = range(1), beta = g.beta, n_points = npts)
    trace(g, tr_g)
    #for s, ls in zip(['up', 'dn'], ['--', ':']):
    #    oplot(g[s], name=i_pp, mode = 'R', ls = ls, color = colors[int(i_pp)])
    plt.plot([x for x in tr_g.mesh], tr_g.data[:,0,0].real, label='$'+i_pp+'$', color = colors[int(i_pp)])
    print i_pp, tr_g.tail[1][0, 0]

g = arch['G_tau']
g = BlockGf(name_block_generator = [(s, rebinning_tau(b, npts)) for s, b in g])
trace(g, tr_g)
plt.plot([x for x in tr_g.mesh], tr_g.data[:,0,0].real, label='$\mathrm{Tr}G_{imp}$', color = 'black')
print 'tot', tr_g.tail[1][0, 0]

fullgbygpp = BlockGf(name_block_generator = [(s, rebinning_tau(b, npts)) for s, b in fullgbygpp])
tr_g = GfImTime(indices = range(1), beta = fullgbygpp.beta, n_points = npts)
trace(fullgbygpp, tr_g)
plt.plot([x for x in tr_g.mesh], tr_g.data[:,0,0].real, label='$sum\;i_{pp}$', ls = '--', color = 'gray')
print "summed pp", tr_g.tail[1][0, 0]

plt.gca().legend(loc = 'lower center', title = "$i_{pp}$")
plt.gca().set_ylim(-1,0)
pp.savefig(plt.gcf(), dpi=300)
pp.close()
