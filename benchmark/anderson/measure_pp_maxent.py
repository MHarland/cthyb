#!/bin/env pytriqs

from maxent.bryanToTRIQS import MaximumEntropy
from pytriqs.archive import HDFArchive
from pytriqs.gf.local import GfReFreq, GfLegendre, GfImTime, rebinning_tau, BlockGf
import numpy as np


par = {"ntau": 100,
       "nomega": 200,
       "bandwidth": 8,
       "sigma": 0.001}
archive_name = "measure_pp.h5"

def trace(g, tr_g):
    tr_g.zero()
    indices = [(b, i, j) for b, i , j in g.all_indices if i == j]
    for ind in indices:
        tr_g += g[ind[0]][ind[1], ind[2]] / len(indices)

def maxent(g, par):
    maxent = MaximumEntropy(g, par['ntau'])
    if par['sigma']:
        maxent.calculateTotDOS(par['nomega'], par['bandwidth'], par['sigma'])
    else:
        maxent.calculateTotDOS(par['nomega'], par['bandwidth'])
    w = maxent.getOmegaMesh()
    a = maxent.getTotDOS()
    return w, a

archive = HDFArchive(archive_name, 'a')
n_pp = len([key for key in archive['G_pp_tau'].keys()])
for i in range(n_pp):
    if i in [0,2,3]: continue
    g = archive['G_pp_tau'][str(i)]
    mesh, a_w = maxent(g, par)
    if i == 1:
        if archive.is_group('maxent'):
            del archive['maxent']
        archive.create_group('maxent')
    archive['maxent'].create_group(str(i))
    res = archive['maxent'][str(i)]
    res['a_w'] = a_w
    res['w'] = mesh

g = archive['G_tau']
mesh, a_w = maxent(g, par)
if archive['maxent'].is_group('tot'):
    del archive['maxent']['tot']
archive['maxent'].create_group('tot')
res = archive['maxent']['tot']
res['a_w'] = a_w
res['w'] = mesh
