from cython cimport typeof

from defs cimport BfSize
from types cimport BfMat, BfType, BfVec

cdef extern from "bf/mat.h":
    BfType bfMatGetType(const BfMat *mat)
    BfSize bfMatGetNumRows(const BfMat *mat)
    BfSize bfMatGetNumCols(const BfMat *mat)
    void bfMatAddInplace(BfMat *mat, const BfMat *otherMat)
    BfMat *bfMatSub(const BfMat *, const BfMat *)
    BfMat *bfMatMul(const BfMat *, const BfMat *)
    BfVec *bfMatMulVec(const BfMat *, const BfVec *)
    BfMat *bfMatRmul(const BfMat *, const BfMat *)
    BfMat *bfMatToType(const BfMat *mat, BfType type)
    void bfMatTranspose(BfMat *mat)
