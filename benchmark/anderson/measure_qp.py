#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *


spin_names = ("up","dn")
mkind = lambda spin: (spin,0) if use_blocks else ("tot",spin)
beta = 10.0
U = 2.0
mu = 1.0
h = 0.1
V = 0.5
epsilon = 2.3
use_blocks = True
use_qn = True
n_iw = 1025
n_tau = 10001
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 5 * 10**4
p["n_cycles"] = 10**6
p["measure_g_qp_tau"] = True
p["use_norm_as_weight"] = True
p["measure_density_matrix"] = True
results_file_name = "measure_qp.h5"

H = U*n(*mkind("up"))*n(*mkind("dn"))
QN = []
if use_qn:
    for spin in spin_names: QN.append(n(*mkind(spin)))
    p["quantum_numbers"] = QN
    p["partition_method"] = "quantum_numbers"
gf_struct = {}
for spin in spin_names:
    bn, i = mkind(spin)
    gf_struct.setdefault(bn,[]).append(i)

mpi.report("Constructing the solver...")
S = Solver(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

mpi.report("Preparing the hybridization function...")
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w << (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon)
for spin in spin_names:
    bn, i = mkind(spin)
    S.G0_iw[bn][i,i] << inverse(iOmega_n + mu - {'up':h,'dn':-h}[spin] - delta_w)

mpi.report("Running the simulation...")
S.solve(h_int=H, **p)

if mpi.is_master_node():
    static_observables = {'Nup' : n(*mkind("up")), 'Ndn' : n(*mkind("dn")), 'unity' : Operator(1.0)}
    with HDFArchive(results_file_name,'w') as Results:
        Results.create_group('G_qp_tau')
        for i, g in enumerate(S.G_qp_tau):
            Results['G_qp_tau'][str(i)] = g
        Results['G_tau'] = S.G_tau
