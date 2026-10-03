#!/usr/bin/env python
'''
Usage: mkofflineveldaily.py <rundir> <firstiter> <enditer> <iterstep>
'''
import sys
import os
import numpy as np
from MITgcmutils import rdmds

args = sys.argv[1:]
rdir = args.pop(0)
it0 = int(args.pop(0))
ite = int(args.pop(0))
itd = int(args.pop(0))

ddir = rdir + '/' + 'diags'
its = range(it0, ite, itd)

os.makedirs(ddir + '/DDuveloff', exist_ok=True)
os.makedirs(ddir + '/DDvveloff', exist_ok=True)

h0fW = rdmds(rdir + '/hFacW')
h0fS = rdmds(rdir + '/hFacS')
mW = h0fW != 0
mS = h0fS != 0

for it in its:
    u = rdmds(ddir + '/DDuvelmass/DDuvelmass', it)
    u[mW] /= h0fW[mW]
    u.astype('>f4').tofile(ddir + '/DDuveloff/DDuveloff.{:010d}.data'.format(it))

    v = rdmds(ddir + '/DDvvelmass/DDvvelmass', it)
    v[mS] /= h0fS[mS]
    v.astype('>f4').tofile(ddir + '/DDvveloff/DDvveloff.{:010d}.data'.format(it))


