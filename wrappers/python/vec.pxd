from defs cimport BfSize
from types cimport BfType, BfVec

cdef extern from "bf/vec.h":
    void bfVecDelete(BfVec **)
    BfType bfVecGetType(BfVec *)
    BfSize bfVecGetSize(const BfVec *)
