#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from mpl_to_latex.matplotlib_to_latex import set_prl_parameters
set_prl_parameters(72, 390)


arch = HDFArchive('measure_pp.h5', 'r')
arch_som = arch['som']
for i_pp, pp in arch_som.items():
    g_w = pp['g_w']
    plt.plot([x for x in g_w.mesh], -1*g_w.data[:,0,0].imag/np.pi, label='$'+i_pp+'$')
plt.gca().set_ylim(bottom = 0)
plt.gca().legend()
plt.savefig("a_pp.pdf", dpi=300)
plt.close()

arch = HDFArchive('measure_pp.h5', 'r')
arch_som = arch['som']
n_graphs = len([k for k in arch_som.keys()])
colors = [plt.cm.jet(i/float(n_graphs)) for i in range(n_graphs)]
for (i_pp, pp), col in zip(arch_som.items(), colors):
    g = pp['g_rec']
    gin = pp['g_in']
    plt.plot([x for x in gin.mesh], gin.data[:,0,0].real, ls = 'dashed', color = col)
    plt.plot([x for x in g.mesh], g.data[:,0,0].real, label='$'+i_pp+'$', color = col)
plt.gca().set_ylim(-1, .05)
plt.gca().legend(fontsize = 10, loc =  'lower center')
plt.savefig('measure_pp_som_g_rec.pdf')
plt.close()
