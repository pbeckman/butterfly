# lbo_MFEM

This directory contains some scripts and executables which use MFEM to set up stiffness and mass matrices for Laplace-Beltrami operator problems.

- `convert_to_uint64_bin.py`: MFEM CSR format uses `int`s for its indices, while butterfly uses `size_t`s. This script uses numpy to convert a binary file containing `int`s to one containing `size_t`s.
- `lbo_MFEM`: Driver program that uses MFEM to set up stiffness and mass matrices.
- `make_plots.py`: Debugging script which loads a stiffness and mass matrix, computes an LBO eigenvector, and plots it.
- `make_sphere_data.sh`: Collects a bunch of h-, p-, and hp-refined stiffness and mass matrices for the sphere.
