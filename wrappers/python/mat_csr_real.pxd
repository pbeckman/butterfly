from size_array cimport BfSizeArray
from types cimport BfMat, BfMatCsrReal, BfTrimesh, BfVec

cdef extern from "bf/mat_csr_real.h":
    BfMatCsrReal *bfMatCsrRealNewViewFactorMatrixFromTrimesh(const BfTrimesh *trimesh, const BfSizeArray *rowInds, const BfSizeArray *colInds)
    BfMat *bfMatCsrRealToMat(BfMatCsrReal *matCsrReal)
