import colorcet as cc
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import scipy.sparse
import sys

from matplotlib.colors import LogNorm
from pathlib import Path

path = Path(sys.argv[1])

L_data = np.fromfile(path/'L_data.bin', np.double)
L_colind = np.fromfile(path/'L_colind.bin', np.uintp)
L_rowptr = np.fromfile(path/'L_rowptr.bin', np.uintp)
L = scipy.sparse.csr_matrix((L_data, L_colind, L_rowptr))

M_data = np.fromfile(path/'M_data.bin', np.double)
M_colind = np.fromfile(path/'M_colind.bin', np.uintp)
M_rowptr = np.fromfile(path/'M_rowptr.bin', np.uintp)
M = scipy.sparse.csr_matrix((M_data, M_colind, M_rowptr))

nodes = np.fromfile(path/'nodes.bin', np.double).reshape(3, -1).T

lam_max = scipy.sparse.linalg.eigsh(L, 1, M, which='LM')[0][0]
print(f'maximum eigenvalue: {lam_max}')

mu = -0.001
Lam, U = scipy.sparse.linalg.eigsh(L, 10, M, mu, which='LM')
u = U[:, 9]

poly_data = pv.PolyData(nodes)
poly_data['u'] = u
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(poly_data, cmap=cc.cm.gouldian, point_size=10)

x, y, z = nodes.T
phi = np.arccos(z)
theta = np.mod(np.arctan2(y, x), 2*np.pi)

plt.figure(figsize=(10, 6))
plt.scatter(theta, phi, 10, u, cmap=cc.cm.gouldian)
plt.xlim(0, 2*np.pi)
plt.ylim(np.pi, 0)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
plt.yticks([np.pi, np.pi/2, 0])
plt.gca().set_xticklabels(['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
plt.gca().set_yticklabels(['0', r'$\pi/2$', r'$\pi$'][::-1])
plt.xlabel(r'$\theta$')
plt.ylabel(r'$\phi$')
plt.show()

# plt.figure(figsize=(11, 10))
# plt.title('L')
# plt.imshow(abs(L.A), cmap=cc.cm.gouldian, norm=LogNorm(clip=True))
# plt.colorbar()
# plt.gca().set_aspect('equal')
# plt.tight_layout()
# plt.show()

# plt.figure(figsize=(11, 10))
# plt.title('M')
# plt.imshow(abs(M.A), cmap=cc.cm.gouldian, norm=LogNorm(clip=True))
# plt.colorbar()
# plt.gca().set_aspect('equal')
# plt.tight_layout()
# plt.show()
