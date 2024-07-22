#pragma once

#include "fac.h"

/* A `FacSpan` is an array of column butterfly factorizations which
 * are concatenated together in a row. This type exists in order to
 * bundle together (possibly) partial butterfly factorizations as
 * they're being streamed column by column.
 *
 * For example, if we have a binary column tree and we're
 * incrementally constructing a butterfly factorization from left to
 * right, we may not want to stream "a power-of-two"-many column
 * blocks, in which case we can't push the butterfly "all the way
 * through". This will leave a remainder of column blocks, and
 * `FacSpan` can be used to group the results together. */
struct BfFacSpan {
  BfSize numFacs;
  BfFac **fac;
};

BfFacSpan *bfFacSpanNewFromPtrArray(BfPtrArray const *ptrArray);
void bfFacSpanInitFromPtrArray(BfFacSpan *facSpan, BfPtrArray const *ptrArray);
void bfFacSpanDeinit(BfFacSpan *facSpan);
void bfFacSpanDealloc(BfFacSpan **facSpan);
void bfFacSpanDelete(BfFacSpan **facSpan);
BfMat *bfFacSpanGetMat(BfFacSpan const *facSpan, BfPolicy policy);
