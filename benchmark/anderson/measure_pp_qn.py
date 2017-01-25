#!/bin/env pytriqs

import numpy as np
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import AtomDiag, act


observables = {'n_tot': n('up', 0) + n('dn', 0), 'n_up': n('up', 0)}

arch = HDFArchive('measure_pp.h5', 'r')
hdiag = arch['h_loc_diagonalization']
dm = arch['density_matrix']
for i_pp in range(hdiag.full_hilbert_space_dim):
    line = 'i_pp: '+str(i_pp)
    state = np.zeros([hdiag.full_hilbert_space_dim])
    state[i_pp] = 1
    for obname, obop in observables.items():
        line += ', '+obname+': '+str(state.dot(act(obop, state, hdiag)))
    print line
