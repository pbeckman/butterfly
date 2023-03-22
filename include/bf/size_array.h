#pragma once

#include "def.h"
#include "types.h"

struct BfSizeArray {
  BfSize *data;
  BfSize size;
  BfSize capacity;
};

BfSizeArray *bfSizeArrayNew();
void bfSizeArrayInitWithDefaultCapacity(BfSizeArray *sizeArray);
void bfSizeArrayDeinit(BfSizeArray *sizeArray);
void bfSizeArrayDealloc(BfSizeArray **sizeArray);
void bfSizeArrayDeinitAndDealloc(BfSizeArray **sizeArray);
void bfSizeArrayExpandCapacity(BfSizeArray *sizeArray, BfSize newCapacity);
void bfSizeArrayAppend(BfSizeArray *sizeArray, BfSize elt);
bool bfSizeArrayContains(BfSizeArray const *sizeArray, BfSize elt);
bool bfSizeArrayIsSorted(BfSizeArray const *sizeArray);
void bfSizeArrayInsertSorted(BfSizeArray *sizeArray, BfSize elt);
BfSize bfSizeArrayFindFirst(BfSizeArray const *sizeArray, BfSize elt);
BfSize bfSizeArrayGet(BfSizeArray const *sizeArray, BfSize i);
void bfSizeArrayCopyData(BfSizeArray const *sizeArray, BfSize *dst);
BfSize bfSizeArrayGetSize(BfSizeArray const *sizeArray);