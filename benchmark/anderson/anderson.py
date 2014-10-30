#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators.operators2 import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

print_master("Welcome to Anderson (1 correlated site + symmetric bath) test.")

H = U*n(*mkind("up"))*n(*mkind("dn"))

QN = []
if use_qn:
    for spin in spin_names: QN.append(n(*mkind(spin)))

p["h_loc"] = H
p["quantum_numbers"] = QN
p["static_observables"] = {'Nup' : n(*mkind("up")), 'Ndn' : n(*mkind("dn"))}

gf_struct = {}
for spin in spin_names:
    bn, i = mkind(spin)
    gf_struct.setdefault(bn,[]).append(i)

print_master("Constructing the solver...")

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

print_master("Preparing the hybridization function...")

# Set hybridization function    
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w <<= (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon)
for spin in spin_names:
    bn, i = mkind(spin)
    S.G0_iw[bn][i,i] <<= inverse(iOmega_n + mu - {'up':h,'dn':-h}[spin] - delta_w)

print_master("Running the simulation...")

# Solve the problem
S.solve(**p)

print "Nup:", S.static_observables["Nup"]
print "Ndn:", S.static_observables["Ndn"]

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name(use_blocks,use_qn),'w')
    for b in gf_struct: Results[b] = S.G_tau[b]
