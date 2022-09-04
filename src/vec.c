#include <bf/vec.h>

#include <assert.h>
#include <math.h>

BfPtr bfVecGetEltPtr(BfVec *x, BfSize i) {
  return (BfByte *)x->data + i*x->stride;
}

BfReal bfVecDist(BfVec const *x, BfVec const *y) {
  assert(x->dtype == y->dtype);

  BfSize n = x->size;
  assert(n == y->size);

  BfReal sum = 0;

  BfByte *x_ptr = x->data, *y_ptr = y->data;

  if (x->dtype == BF_DTYPE_REAL) {
    BfReal diff = 0;
    for (BfSize i = 0; i < n; ++i) {
      diff = *(BfReal *)x_ptr - *(BfReal *)y_ptr;
      sum += diff*diff;
      x_ptr += x->stride;
      y_ptr += y->stride;
    }
  }

  if (x->dtype == BF_DTYPE_COMPLEX) {
    BfComplex diff = 0;
    for (BfSize i = 0; i < n; ++i) {
      diff = *(BfComplex *)x_ptr - *(BfComplex *)y_ptr;
      sum += diff*conj(diff);
      x_ptr += x->stride;
      y_ptr += y->stride;
    }
  }

  return sqrt(sum);
}

void bfVecScaleByReal(BfVec *x, BfReal scale) {
  BfSize n = x->size;
  BfByte *x_ptr = x->data;

  if (x->dtype == BF_DTYPE_REAL) {
    for (BfSize i = 0; i < n; ++i) {
      *(BfReal *)x_ptr *= scale;
      x_ptr += x->stride;
    }
  }

  if (x->dtype == BF_DTYPE_COMPLEX) {
    for (BfSize i = 0; i < n; ++i) {
      *(BfComplex *)x_ptr *= scale;
      x_ptr += x->stride;
    }
  }
}