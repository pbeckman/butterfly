#include <bf/mat_csr_real.h>
#include <bf/size_array.h>
#include <bf/trimesh.h>
#include <bf/util.h>

#ifndef BF_EMBREE
#  error "Can only build examples/radiosity when compiled with Embree"
#endif

#include <stdlib.h>

int main(void) {
  bfToc();

  char const *path = "67p.obj";

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(path);

  if (!bfTrimeshHasVertexNormals(trimesh)) {
    printf("tried to load obj file (\"%s\") without vertex normals", path);
    exit(EXIT_FAILURE);
  }

  if (bfTrimeshHasVertexNormals(trimesh) && !bfTrimeshHasFaceNormals(trimesh))
    bfTrimeshComputeFaceNormalsMatchingVertexNormals(trimesh);

  bfTrimeshInitEmbree(trimesh);

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  BfSize numFaces = bfTrimeshGetNumFaces(trimesh);
  printf("loaded mesh with %lu verts and %lu faces [%0.2fs]\n", numVerts, numFaces, bfToc());

  BfSizeArray *rowInds = bfSizeArrayNewIota(numFaces);
  BfSizeArray *colInds = bfSizeArrayNewIota(numFaces);
  BfMatCsrReal *mat = bfMatCsrRealNewViewFactorMatrixFromTrimesh(trimesh, rowInds, colInds);

  printf("computed view factor matrix [%0.2fs]\n", bfToc());

  bfMatCsrRealDeinitAndDealloc(&mat);

  return EXIT_SUCCESS;
}
