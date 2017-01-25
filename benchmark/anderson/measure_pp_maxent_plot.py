#!/bin/env pytriqs

from pytriqs.archive import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt


pp = PdfPages('a_pp_maxent.pdf')
plt.clf()
arch = HDFArchive('measure_pp.h5', 'r')
arch_m = arch['maxent']
for i_pp, pp in arch_m.items():
    a_w = pp['a_w']
    mesh = pp['w']
    plt.plot(mesh, a_w, label=str(i_pp))
a_w = arch_m['tot']['a_w']
mesh = pp['w']
plt.plot(mesh, a_w, label='tot')
plt.gca().set_ylim(bottom = 0)
pp.savefig(plt.gcf())
pp.close()
