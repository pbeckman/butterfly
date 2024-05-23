#!/usr/bin/env python

import colorcet as cc
import numpy as np
import pyvista as pv
import sys

pv.set_plot_theme('document')

depth = int(sys.argv[1])

grid = pv.read(sys.argv[2])
perm = np.fromfile('perm.bin', dtype=np.uintp)

verts = grid.points.copy()
faces = grid.faces.reshape(-1, 4)[:, 1:]

n = perm.size
assert n == verts.shape[0]

rev_perm = np.empty_like(perm)
for i in range(n):
    rev_perm[perm[i]] = i

rev_perm = n*np.floor(8*rev_perm/n)/8

with open('node_spans.txt', 'r') as f:
    D, I0, I1 = \
        np.array([[int(_) for _ in line.split()] for line in f.readlines()]).T

scalars = np.empty_like(rev_perm)
for k, (i0, i1) in enumerate(zip(I0[D == depth], I1[D == depth])):
    scalars[perm[i0:i1]] = k

poly_data = pv.PolyData(verts)
poly_data['scalars'] = scalars

_ = pv.Plotter()
_.add_mesh(grid, color='white')
_.add_mesh(
    poly_data,
    render_points_as_spheres=True,
    cmap=cc.cm.glasbey,
    # cmap='cool',
    smooth_shading=False,
    show_scalar_bar=False,
    point_size=9)
_.show()
