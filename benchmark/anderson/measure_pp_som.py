#!/bin/env pytriqs

from pytriqs.archive import HDFArchive
from pytriqs.gf.local import GfReFreq, GfLegendre, GfImTime, rebinning_tau, BlockGf
from pytriqs.utility import mpi
from triqs_som.som import Som
import numpy as np


npts = 100
run_params = {}
run_params['energy_window'] = (-3, 3)
run_params['max_time'] = -1
run_params['verbosity'] = 2
run_params['t'] = 1000
run_params['f'] = 50
run_params['adjust_f'] = False
run_params['l'] = 100
run_params['adjust_l'] = False
run_params['make_histograms'] = False
run_params['hist_max'] = 2
run_params['hist_n_bins'] = 100
archive_name = "measure_pp.h5"

def trace(g, tr_g):
    tr_g.zero()
    indices = [(b, i, j) for b, i , j in g.all_indices if i == j]
    for ind in indices:
        tr_g += g[ind[0]][ind[1], ind[2]] / len(indices)

def som(g, npts, run_params):
    if npts is not None: g = BlockGf(name_block_generator = [(s, rebinning_tau(b, npts)) for s, b in g])
    npts = len([x for x in g.mesh])
    tr_g = GfImTime(indices = range(1), beta = g.beta, n_points = npts)
    trace(g, tr_g)
    s = tr_g.copy()
    g_rec = tr_g.copy()
    gw = GfReFreq(window = (run_params['energy_window'][0], run_params['energy_window'][1]), n_points = 5000, indices = tr_g.indices)
    som = Som(tr_g, s, "FermionGf", tr_g.tail[1][0, 0].real)
    som.run(**run_params)
    g_rec << som
    gw << som
    return tr_g, g_rec, gw, s

n_pp = None
if mpi.is_master_node():
    archive = HDFArchive(archive_name, 'a')
    n_pp = len([key for key in archive['G_pp_tau'].keys()])
n_pp = mpi.bcast(n_pp)

for i in range(n_pp):
    g = None
    if mpi.is_master_node():
        g = archive['G_pp_tau'][str(i)]
    g = mpi.bcast(g)
    tr_g, g_rec, gw, s = som(g, npts, run_params)
    if mpi.is_master_node():
        if i == 0:
            if archive.is_group('som'):
                del archive['som']
            archive.create_group("som")
        archive["som"].create_group(str(i))
        res = archive["som"][str(i)]
        res['g_in'] = tr_g
        res['g_rec'] = g_rec
        res['g_w'] = gw
        res['s'] = s
        #res['histograms'] = som.histograms
        res['parameters'] = run_params

g = None
if mpi.is_master_node():
    g = archive['G_tau']
g = mpi.bcast(g)
tr_g, g_rec, gw, s = som(g, npts, run_params)
if mpi.is_master_node():
    if archive['som'].is_group('tot'):
        del archive['som']['tot']
    archive['som'].create_group('tot')
    res = archive['som']['tot']
    res['g_in'] = tr_g
    res['g_rec'] = g_rec
    res['g_w'] = gw
    res['s'] = s
    #res['histograms'] = som.histograms
    res['parameters'] = run_params
