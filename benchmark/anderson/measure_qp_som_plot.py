#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


pp = PdfPages('a_qp.pdf')
plt.clf()
arch = HDFArchive('measure_qp.h5', 'r')
arch_som = arch['som']
for i_qp, qp in arch_som.items():
    g_w = qp['g_w']
    oplot(g_w, name=i_qp, mode = 'S')
g_w = arch_som['tot']['g_w']
oplot(g_w * .5, name=i_qp, mode = 'S')
plt.gca().set_ylim(bottom = 0)
pp.savefig(plt.gcf())
pp.close()
