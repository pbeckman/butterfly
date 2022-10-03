#include <bf/mat_sum.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatSum)
#undef INTERFACE

BF_STUB(BfMat *, MatSumCopy, BfMat const *)
BF_STUB(BfMat *, MatSumGetView, BfMat *)
BF_STUB(BfVec *, MatSumGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatSumGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatSumGetColRangeView, BfMat *, BfSize, BfSize, BfSize)
BF_STUB(void, MatSumDelete, BfMat **)
BF_STUB(BfMat *, MatSumEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatSumZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatSumGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_SUM;
}

BF_STUB(BfSize, MatSumNumBytes, BfMat const *)
BF_STUB(void, MatSumSave, BfMat const *, char const *)
BF_STUB(void, MatSumPrint, BfMat const *, FILE *)

BfSize bfMatSumGetNumRows(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = BF_SIZE_BAD_VALUE;

  BfMatSum const *matSum = bfMatConstToMatSumConst(mat);
  HANDLE_ERROR();

  numRows = bfMatGetNumRows(bfMatSumGetTermConst(matSum, 0));

  END_ERROR_HANDLING() {}

  return numRows;
}

BfSize bfMatSumGetNumCols(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfSize numCols = BF_SIZE_BAD_VALUE;

  BfMatSum const *matSum = bfMatConstToMatSumConst(mat);
  HANDLE_ERROR();

  numCols = bfMatGetNumCols(bfMatSumGetTermConst(matSum, 0));

  END_ERROR_HANDLING() {}

  return numCols;
}

BF_STUB(void, MatSumSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatSumSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatSumSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatSumGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatSumGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatSumGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatSumGetColRangeCopy , BfMat const *, BfSize, BfSize)
BF_STUB(void, MatSumSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatSumPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatSumPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatSumRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatSumColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatSumColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatSumColNorms, BfMat const *)
BF_STUB(void, MatSumScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatSumSumCols, BfMat const *)
BF_STUB(void, MatSumAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatSumAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatSumSub, BfMat const *, BfMat const *)
BF_STUB(void, MatSumSubInplace, BfMat *, BfMat const *)

BfMat *bfMatSumMul(BfMat const *mat, BfMat const *otherMat) {
  BEGIN_ERROR_HANDLING();

  BfMatSum const *matSum = NULL;
  BfMat *result = NULL;

  BfSize m = bfMatGetNumRows(mat);
  BfSize p = bfMatGetNumCols(mat);
  BfSize n = bfMatGetNumCols(otherMat);

  matSum = bfMatConstToMatSumConst(mat);
  HANDLE_ERROR();

  if (p != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  result = bfMatZerosLike(otherMat, m, n);

  for (BfSize i = 0; i < bfMatSumGetNumTerms(matSum); ++i) {
    BfMat const *term = bfMatSumGetTermConst(matSum, i);
    BfMat *prod = bfMatMul(term, otherMat);
    bfMatAddInplace(result, prod);
    bfMatDelete(&prod);
  }

  END_ERROR_HANDLING()
    bfMatDelete(&result);

  return result;
}

BF_STUB(BfVec *, MatSumMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatSumMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatSumSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatSumLstSq, BfMat const *, BfMat const *)
BF_STUB(bool, MatSumIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatSumBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatSumIsZero, BfMat const *)

/** Upcasting: */

BfMat *bfMatSumToMat(BfMatSum *matSum) {
  return &matSum->super;
}

/** Downcasting: */

BfMatSum const *bfMatConstToMatSumConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_SUM)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatSum const *)mat;
  }
}

/** Implementation: MatSum */

BfMatSum *bfMatSumNew() {
  BEGIN_ERROR_HANDLING();

  BfMatSum *sum = malloc(sizeof(BfMatSum));
  if (sum == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return sum;
}

void bfMatSumInit(BfMatSum *mat) {
  BEGIN_ERROR_HANDLING();

  /* We don't store the number of rows or columns in `mat->super`
   * since we always look up the number of rows and columns from the
   * terms at runtime. */
  bfMatInit(&mat->super, &MatVtbl, BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE);

  mat->termArr = bfGetUninitializedPtrArray();

  bfInitPtrArray(&mat->termArr, /* capacity: */ 4);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfPtrArrayDeinit(&mat->termArr);
}

void bfMatSumDeinit(BfMatSum *sum) {
  for (BfSize i = 0; i < bfPtrArraySize(&sum->termArr); ++i) {
    BfMat *mat = bfPtrArrayGet(&sum->termArr, i);
    bfMatDelete(&mat);
  }

  bfPtrArrayDeinit(&sum->termArr);
}

void bfMatSumDealloc(BfMatSum **sum) {
  free(*sum);
  *sum = NULL;
}

void bfMatSumDeinitAndDealloc(BfMatSum **sum) {
  bfMatSumDeinit(*sum);
  bfMatSumDealloc(sum);
}

BfSize bfMatSumGetNumTerms(BfMatSum const *sum) {
  return bfPtrArraySize(&sum->termArr);
}

void bfMatSumAddTerm(BfMatSum *sum, BfMat *term) {
  bfPtrArrayAppend(&sum->termArr, term);
}

BfMat *bfMatSumGetTerm(BfMatSum *sum, BfSize i) {
  return bfPtrArrayGet(&sum->termArr, i);
}

BfMat const *bfMatSumGetTermConst(BfMatSum const *sum, BfSize i) {
  return bfPtrArrayGet(&sum->termArr, i);
}