#!/usr/bin/env python3

import sys; sys.path.insert(-1, '../../wrappers/python')

import argparse
import matplotlib.pyplot as plt
import numpy as np

import butterfly as bf

# simulation parameters
xmin, xmax = -1, 1
ymin, ymax = -1, 1
r = 0.2 # spacing between ellipses
a, b = 0.04, 0.08 # range of ellipse semi axes
h = 0.02 # point spacing

# seed bf's PRNG
bf.seed(0)

# use Poisson disk sampling to generate ellipse centers
bbox = bf.Bbox2(xmin, xmax, ymin, ymax)
centers = bf.Points2.sample_poisson_disk(bbox, r)
num_ellipses = len(centers)

# randomly sample ellipses
ellipses = []
for i in range(num_ellipses):
    a, b = np.random.uniform((a, b), size=2)
    theta = np.random.uniform(2*np.pi)
    ellipses.append(bf.Ellipse(max(a, b), min(a, b), centers[i], theta))

# discretize the ellipses to get the problem geometry
X, W, N = bf.Points2(), bf.RealArray(), bf.Vectors2()
for ellipse in ellipses:
    n = int(np.ceil(ellipse.perimeter/h))
    X_, _, N_, W_ = ellipse.sample_linspaced(n)
    X.extend(X_)
    N.extend(N_)
    W.extend(W_)

# build quadtree on points w/ normals
quadtree = bf.Quadtree(X, N)