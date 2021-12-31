#pragma once

#include "def.h"
#include "error.h"

static BfSize const BF_ARRAY_DEFAULT_CAPACITY = 1024;

enum BfPtrArrayFlags {
  BF_PTR_ARRAY_FLAG_NONE = 0,
  BF_PTR_ARRAY_FLAG_VIEW = 1 << 0
};

typedef struct BfPtrArray {
  enum BfPtrArrayFlags flags;
  BfPtr *data;
  BfSize capacity, num_elts;
} BfPtrArray;

typedef void (*BfPtrFunc)(BfPtr elt_ptr, BfPtr arg_ptr);

BfPtrArray bfGetUninitializedPtrArray();

void
bfInitPtrArray(BfPtrArray *arr, BfSize capacity);

void
bfInitPtrArrayWithDefaultCapacity(BfPtrArray *arr);

void bfMakeEmptyPtrArrayView(BfPtrArray *arr);

void
bfFreePtrArray(BfPtrArray *arr);

BfSize bfPtrArraySize(BfPtrArray const *arr);

bool bfPtrArrayIsEmpty(BfPtrArray const *arr);

void
bfPtrArrayAppend(BfPtrArray *arr, BfPtr ptr);

void
bfPtrArrayGet(BfPtrArray const *arr, BfSize pos, BfPtr *ptr);

void
bfPtrArrayGetFirst(BfPtrArray const *arr, BfPtr *ptr);

void
bfPtrArrayGetLast(BfPtrArray const *arr, BfPtr *ptr);

void
bfPtrArrayGetRangeView(BfPtrArray const *arr, BfSize start, BfSize end,
                       BfPtrArray *view);

void
bfMapPtrArray(BfPtrArray *arr, BfPtrFunc func, BfPtr ptr);
