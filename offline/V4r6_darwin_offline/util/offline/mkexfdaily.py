#!/usr/bin/env python
'''
Usage: mkexfdaily.py <diagdir> <years>...

Ex: ./mkexfdaily.py runonline/diags/DSIarea {1992..2025}

compiles runonline/diags/DSIarea/DSIarea.########## into yearly files
for the given years, written to runonline/diags/DSIarea_####.
'''
import sys
from os.path import split as psplit, join as pjoin
import numpy as np
from datetime import datetime, timedelta
from MITgcmutils import rdmds, wrmds

# from data.cal
startdate = datetime(1992, 1, 1, 12)
# from data
niter0 = 1
deltaT = 3600.
# averaging frequency
freq = 86400.
# averages are labelled at end of the interval, need to be placed at center
labeloffset = freq/2

args = sys.argv[1:]
ddir = args.pop(0)

_, fld = psplit(ddir)
prefix = pjoin(ddir, fld)
#prefix = pjoin(ddir, '_')

for arg in args:
    year = int(arg)
    d0 = datetime(year, 1, 1)
    de = datetime(year+1, 1, 1)
    t0 = (d0 - startdate).total_seconds()
    te = (de - startdate).total_seconds()
    ts = np.arange(t0, te, freq) + labeloffset
    its = np.round(ts/deltaT).astype(int)

    print(startdate + its[0]*timedelta(seconds=deltaT), its[0],
          startdate + its[-1]*timedelta(seconds=deltaT), its[-1])

    if t0 < niter0*freq:
        # missing an average before run starts
        ain = rdmds(prefix, list(its[1:]), astype='>f4')
        a = np.zeros((len(its),) + ain.shape[1:], '>f4')
        # copy first available average
        a[0] = ain[0]
        a[1:] = ain
    else:
        a = rdmds(prefix, list(its), astype='>f4')

    a.tofile(ddir+'_{}'.format(year))

