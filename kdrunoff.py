import sys
import numpy as np
from scipy import spatial

def remap_runoff(runoff, land_sea_mask, x_t, y_t, spread=1):
    """
    Return a runoff field with all values moved from land to the nearest ocean
    point.

    The land sea mask is 0 on land points and 1 elsewhere.

    Instead of redistributing to a single nearest point the runoff can be
    spread to a number.
    """

    output = np.copy(runoff)

    coords = [list(i) for i in zip(x_t.ravel(), y_t.ravel())]

    xv, yv = np.meshgrid(range(x_t.shape[1]), range(x_t.shape[0]))
    indices = [list(i) for i in zip(xv.ravel(), yv.ravel())]

    land_points = []
    land_indices = []
    ocean_points = []
    ocean_indices = []
    for n, i in enumerate(indices):
        if land_sea_mask[i[1], i[0]] < 0.5:
            land_points.append(coords[n])
            land_indices.append(indices[n])
        else:
            ocean_points.append(coords[n])
            ocean_indices.append(indices[n])

    tree = spatial.KDTree(ocean_points)

    # Move runoff from land to ocean.
    for n, i in enumerate(land_points):
        val = runoff[land_indices[n][1], land_indices[n][0]]
        if val > 0.0:
            output[land_indices[n][1], land_indices[n][0]] -= val

            distance, index = tree.query(i)
            output[ocean_indices[index][1], ocean_indices[index][0]] += val

    return output
