#!/bin/env pytriqs

import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.parameters.parameters import Parameters
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
import itertools
from collections import OrderedDict

from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

def print_master(msg):
    if mpi.rank==0: print msg

print_master("Welcome to asymm_bath test (1 band with a small asymmetric hybridization function).")
print_master("This test helps to detect sampling problems.")

# H_loc parameters
beta = 40.0
ed = -1.0
U = 2.0

#epsilon = [0.0,0.1,0.2,0.3,0.4,0.5]
V = 0.2
epsilon = [0.3]

# Parameters
p = Parameters()
p["beta"] = beta
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 20000
p["n_cycles"] = 2000000
p["n_tau_delta"] = 1000
p["n_tau_g"] = 1000
p["measure_pert_order"] = True

# Block structure of GF
gf_struct = OrderedDict()
gf_struct['up'] = [0]
gf_struct['down'] = [0]

# Hamiltonian
H = ed*(N("up",0) + N("down",0)) + U*N("up",0)*N("down",0)

# Quantum numbers
qn = []
qn.append(N("up",0))
qn.append(N("down",0))
    
S = Solver(parameters=p, H_local=H, quantum_numbers=qn, gf_struct=gf_struct)

if mpi.rank==0:
    arch = HDFArchive('asymm_bath.h5','w')
    G_pp = PdfPages('G_asymm_bath.pdf')
    Delta_pp = PdfPages('Delta_asymm_bath.pdf')

# Set hybridization function
for e in epsilon:
    delta_w = GfImFreq(indices = [0], beta=beta)
    delta_w <<= (V**2) * inverse(iOmega_n - e)

    S.Delta_tau["up"] <<= InverseFourier(delta_w)
    S.Delta_tau["down"] <<= InverseFourier(delta_w)

    S.solve(parameters=p)
  
    if mpi.rank==0:
        arch['epsilon_' + str(e)] = {"G_up":S.G_tau["up"], "G_down":S.G_tau["down"],
                                     "Delta_up":S.Delta_tau["up"], "Delta_down":S.Delta_tau["down"]}

        # Plot G(\tau)
        plt.clf()
        oplot(S.G_tau["up"], name="$\uparrow\uparrow$")
        oplot(S.G_tau["down"],name="$\downarrow\downarrow$")
        
        axes = plt.gca()
        axes.set_ylabel('$G(\\tau)$')
        axes.legend(loc='lower center',prop={'size':10})
        
        axes.set_title("$U=%.1f$, $\epsilon_d=%.1f$, $V=%.1f$, $\epsilon_k=%.1f$" % (U,ed,V,e))

        G_pp.savefig(plt.gcf())

        # Plot Delta(\tau)
        plt.clf()
        oplot(S.Delta_tau["up"], name="$\uparrow\uparrow$")
        oplot(S.Delta_tau["down"],name="$\downarrow\downarrow$")

        axes = plt.gca()
        axes.set_ylabel('$\Delta(\\tau)$')
        axes.legend(loc='lower center',prop={'size':10})

        axes.set_title("$V=%.1f$, $\epsilon_k=%.1f$" % (V,e))

        Delta_pp.savefig(plt.gcf())

if mpi.rank==0:
    G_pp.close()
    Delta_pp.close()
