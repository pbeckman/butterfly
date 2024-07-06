#!/usr/bin/env sh

mkdir -p sphere_meshes/h
mkdir -p sphere_meshes/p
mkdir -p sphere_meshes/hp

move_to_dir () {
	mv L_data.bin $1
	mv L_colind.bin $1
	mv L_rowptr.bin $1
	mv M_data.bin $1
	mv M_colind.bin $1
	mv M_rowptr.bin $1
	mv nodes.bin $1
}

convert_single_binary_file_to_uint64 () {
	python <<EOF
import numpy as np
import sys
np.fromfile('${1}', dtype=np.int32).astype(np.uint64).tofile('${1}')
EOF
}

convert_binary_files_to_uint64 () {
	convert_single_binary_file_to_uint64 ${1}/L_colind.bin
	convert_single_binary_file_to_uint64 ${1}/L_rowptr.bin
	convert_single_binary_file_to_uint64 ${1}/M_colind.bin
	convert_single_binary_file_to_uint64 ${1}/M_rowptr.bin
}

# fix h to some coarse level (e.g. ~5k triangles) and take p from 1 up
# to ~20 (or whatever the largest p is that gives ~1M degrees of
# freedom)

for i in $(seq 1 8); do
	./lbo_MFEM --elem 0 --order ${i} --refine 5
	mkdir -p sphere_meshes/p/${i}
	move_to_dir sphere_meshes/p/${i}
	convert_binary_files_to_uint64 sphere_meshes/p/${i}
done

# fix p=1 and for h such that the mesh has ~5k up to ~1M triangles
# (you can replace the existing sphere meshes)

for i in $(seq 4 8); do
	./lbo_MFEM --elem 0 --order 1 --refine ${i}
	mkdir -p sphere_meshes/h/${i}
	move_to_dir sphere_meshes/h/${i}
	convert_binary_files_to_uint64 sphere_meshes/h/${i}
done

# pick moderate parameters like p=6 and h such that there are ~100k
# total DOFs

./lbo_MFEM --elem 0 --order 5 --refine 4
mkdir -p sphere_meshes/hp
move_to_dir sphere_meshes/hp
convert_binary_files_to_uint64 sphere_meshes/hp
