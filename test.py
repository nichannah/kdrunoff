#!/usr/bin/env python

import sys
import numpy as np
import netCDF4 as nc
from kdrunoff import remap_runoff

def main():

    with nc.Dataset('test.nc') as f:
        x_t = f.variables['x_T'][:]
        y_t = f.variables['y_T'][:]
        mask = f.variables['wet'][:]

    runoff = np.random.random(mask.shape)

    # Make sure there is something on the land.
    assert np.sum((1 - mask[:])*runoff[:]) > 0.0

    new_runoff = remap_runoff(runoff, mask, x_t, y_t)

    # Make sure there is nothing left of the land.
    assert np.sum((1 - mask[:])*new_runoff) == 0.0

    # The totals should be the same
    assert np.sum(runoff) == np.sum(new_runoff)

if __name__ == '__main__':
    sys.exit(main())
