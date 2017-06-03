import sys
import numpy as np
from scipy import spatial

class Kdrunoff:

    def __init__(self, land_sea_mask, x_t, y_t):

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

        self.tree = spatial.KDTree(ocean_points)
        self.land_points = land_points
        self.land_indices = land_indices
        self.ocean_indices = ocean_indices

    def remap(self, runoff, spread=1):
        """
        Return a runoff field with all values moved from land to the nearest ocean
        point.

        The land sea mask is 0 on land points and 1 elsewhere.

        Instead of redistributing to a single nearest point the runoff can be
        spread to a number.
        """

        output = np.copy(runoff)

        # Move runoff from land to ocean.
        for n, i in enumerate(self.land_points):
            val = runoff[self.land_indices[n][1], self.land_indices[n][0]]
            if val > 0.0:
                output[self.land_indices[n][1], self.land_indices[n][0]] -= val

                distance, index = self.tree.query(i)
                output[self.ocean_indices[index][1], self.ocean_indices[index][0]] += val

        return output
