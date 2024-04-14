#!/usr/bin/env python

import colorcet as cc
import numpy as np
import pyvista as pv

verts = np.fromfile('verts.bin', dtype=np.float64).reshape(-1, 3)
faces = np.fromfile('faces.bin', dtype=np.uintp).reshape(-1, 3)
phi = np.fromfile('phi.bin', dtype=np.float64)

assert verts.shape[0] == phi.size

abs_phi_max = abs(phi).max()
clim = (-abs_phi_max, abs_phi_max)

grid = pv.make_tri_mesh(verts, faces)
grid['phi'] = phi

_ = pv.Plotter()
_.add_mesh(grid, show_edges=True, cmap=cc.cm.coolwarm, clim=clim)
_.show()
