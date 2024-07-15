from defs cimport BfReal, BfSize
from types cimport BfVec, BfVecReal

cdef extern from "bf/vec_real.h":
    BfVec *bfVecRealToVec(BfVecReal *vecReal)
    BfVecReal *bfVecToVecReal(BfVec *vec)
    BfVecReal *bfVecRealNewViewFromPtr(BfSize size, BfReal *data, BfSize stride)
    BfReal bfVecRealGetElt(const BfVecReal *vecReal, BfSize i)
    BfReal *bfVecRealGetDataPtr(BfVecReal *vecReal)
