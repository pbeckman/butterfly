#include <bf/mat_block_diag.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/indexed_mat.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_zero.h>
#include <bf/mem.h>
#include <bf/ptr_array.h>
#include <bf/util.h>
#include <bf/vec_real.h>
#include <bf/vec_zero.h>

#include "macros.h"

/** Helper macros: */

#define NUM_ROW_BLOCKS(mat) mat->super.super.numRows
#define NUM_COL_BLOCKS(mat) mat->super.super.numCols
#define ROW_OFFSET(mat, k) mat->super.rowOffset[k]
#define COL_OFFSET(mat, k) mat->super.colOffset[k]
#define NUM_BLOCK_ROWS(mat, k) ROW_OFFSET(mat, k + 1) - ROW_OFFSET(mat, k)
#define NUM_BLOCK_COLS(mat, k) COL_OFFSET(mat, k + 1) - COL_OFFSET(mat, k)
#define BLOCK(mat, k) mat->super.block[k]

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatBlockDiagGetView,
  .Copy = (__typeof__(&bfMatBlockDiagCopy))bfMatBlockDiagCopy,
  .Steal = (__typeof__(&bfMatSteal))bfMatBlockDiagSteal,
  .GetRowCopy = (__typeof__(&bfMatBlockDiagGetRowCopy))bfMatBlockDiagGetRowCopy,
  .Delete = (__typeof__(&bfMatBlockDiagDelete))bfMatBlockDiagDelete,
  .GetType = (__typeof__(&bfMatBlockDiagGetType))bfMatBlockDiagGetType,
  .NumBytes = (__typeof__(&bfMatNumBytes))bfMatBlockDiagNumBytes,
  .Dump = (__typeof__(&bfMatDump))bfMatBlockDiagDump,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatBlockDiagGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatBlockDiagGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatBlockDiagGetRowRangeCopy,
  .ScaleCols = (__typeof__(&bfMatBlockDiagScaleCols))bfMatBlockDiagScaleCols,
  .Mul = (__typeof__(&bfMatBlockDiagMul))bfMatBlockDiagMul,
  .MulVec = (__typeof__(&bfMatMulVec))bfMatBlockDiagMulVec,
  .RmulVec = (__typeof__(&bfMatMulVec))bfMatBlockDiagRmulVec,
  .Negate = (__typeof__(&bfMatBlockDiagNegate))bfMatBlockDiagNegate,
  .PrintBlocksDeep = (__typeof__(&bfMatPrintBlocksDeep))bfMatBlockDiagPrintBlocksDeep,
  .Solve = (__typeof__(&bfMatSolve))bfMatBlockDiagSolve,
  .Transpose = (__typeof__(&bfMatTranspose))bfMatBlockDiagTranspose,
};

BfMat *bfMatBlockDiagGetView(BfMatBlockDiag *matBlockDiag) {
  BF_ERROR_BEGIN();

  BfMatBlockDiag *matBlockDiagView = bfMatBlockDiagNew();
  HANDLE_ERROR();

  *matBlockDiagView = *matBlockDiag;

  BfMat *matView = bfMatBlockDiagToMat(matBlockDiagView);

  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return matView;
}

BfMat *bfMatBlockDiagCopy(BfMat const *mat) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = NULL;
  BfMatBlockDiag const *matBlockDiag = NULL;
  BfMatBlockDiag *matBlockDiagCopy = NULL;
  BfMat const *block = NULL;
  BfMat *blockCopy = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockDiag = bfMatConstToMatBlockDiagConst(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);

  matBlockDiagCopy = bfMatBlockDiagNew();
  HANDLE_ERROR();

  bfMatBlockDiagInit(matBlockDiagCopy, numBlocks, numBlocks);
  HANDLE_ERROR();

  for (BfSize i = 0; i <= numBlocks; ++i)
    matBlockDiagCopy->super.rowOffset[i] = matBlock->rowOffset[i];

  for (BfSize j = 0; j <= numBlocks; ++j)
    matBlockDiagCopy->super.colOffset[j] = matBlock->colOffset[j];

  for (BfSize i = 0; i < numBlocks; ++i) {
    block = bfMatBlockDiagGetBlockConst(matBlockDiag, i);
    blockCopy = bfMatCopy(block);
    HANDLE_ERROR();
    bfMatBlockDiagSetBlock(matBlockDiagCopy, i, blockCopy);
  }

  BF_ERROR_END() {
    bfMatDelete(&blockCopy);
    bfMatBlockDiagDeinitAndDealloc(&matBlockDiagCopy);
  }

  return bfMatBlockDiagToMat(matBlockDiagCopy);
}

BfMat *bfMatBlockDiagSteal(BfMatBlockDiag *matBlockDiag) {
  BF_ERROR_BEGIN();

  BfMat *mat = bfMatBlockDiagToMat(matBlockDiag);

  if (bfMatIsView(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockDiag *matBlockDiagNew = bfMatBlockDiagNew();
  HANDLE_ERROR();

  *matBlockDiagNew = *matBlockDiag;

  mat->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return bfMatBlockDiagToMat(matBlockDiagNew);
}

BfVec *bfMatBlockDiagGetRowCopy(BfMat const *mat, BfSize i) {
  BF_ERROR_BEGIN();

  BfVec *rowCopy = NULL;
  BfMatBlock const *matBlock = NULL;
  BfMatBlockDiag const *matBlockDiag = NULL;
  BfMat const *block = NULL;

  matBlock = bfMatConstToMatBlockConst(mat);
  HANDLE_ERROR();

  matBlockDiag = bfMatConstToMatBlockDiagConst(mat);
  HANDLE_ERROR();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (numRowBlocks == 0)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);


  BfSize k = 0, i0 = BF_SIZE_BAD_VALUE, i1 = BF_SIZE_BAD_VALUE;
  for (; k < numRowBlocks; ++k) {
    i0 = matBlockDiag->super.rowOffset[k];
    i1 = matBlockDiag->super.rowOffset[k + 1];
    if (i0 <= i && i < i1)
      break;
  }
  BF_ASSERT(BF_SIZE_OK(i0));
  BF_ASSERT(BF_SIZE_OK(i1));

  block = BLOCK(matBlockDiag, k);

  rowCopy = bfMatGetRowCopy(block, i - i0);
  HANDLE_ERROR();

  /* Pad what we have so far on either side to make sure it matches
   * the width of `mat` */

  BfSize j0 = matBlockDiag->super.colOffset[k];

  BfVec *cat = NULL;
  BfVecZero *zeroPad = NULL;

  if (j0 > 0) {
    zeroPad = bfVecZeroNew();
    HANDLE_ERROR();

    bfVecZeroInit(zeroPad, j0);
    HANDLE_ERROR();

    cat = bfVecConcat(bfVecZeroToVec(zeroPad), rowCopy);
    HANDLE_ERROR();

    bfVecZeroDeinitAndDealloc(&zeroPad);
    bfVecDelete(&rowCopy);

    rowCopy = cat;
  }

  BfSize numCols = bfMatGetNumCols(mat);

  if (rowCopy->size < numCols) {
    zeroPad = bfVecZeroNew();
    HANDLE_ERROR();

    bfVecZeroInit(zeroPad, numCols - rowCopy->size);
    HANDLE_ERROR();

    cat = bfVecConcat(rowCopy, bfVecZeroToVec(zeroPad));
    HANDLE_ERROR();

    bfVecZeroDeinitAndDealloc(&zeroPad);
    bfVecDelete(&rowCopy);

    rowCopy = cat;
  }

  BF_ASSERT(rowCopy->size == numCols);

  BF_ERROR_END() {}

  return rowCopy;
}

void bfMatBlockDiagDelete(BfMat **mat) {
  bfMatBlockDiagDeinitAndDealloc((BfMatBlockDiag **)mat);
}

BfType bfMatBlockDiagGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_BLOCK_DIAG;
}

BfSize bfMatBlockDiagNumBytes(BfMatBlockDiag const *matBlockDiag) {
  BfSize numBytes = 0;
  for (BfSize i = 0; i < bfMatBlockDiagNumBlocks(matBlockDiag); ++i)
    numBytes += bfMatNumBytes(BLOCK(matBlockDiag, i));
  return numBytes;
}

void bfMatBlockDiagDump(BfMatBlockDiag const *matBlockDiag, FILE *fp) {
  BF_ERROR_BEGIN();

  /* Serialize the number of blocks: */
  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);
  fwrite(&numBlocks, sizeof(BfSize), 1, fp);

  // TODO: lift this restriction...
  BF_ASSERT(NUM_ROW_BLOCKS(matBlockDiag) == NUM_COL_BLOCKS(matBlockDiag));

  /* Serialize row and column offsets: */
  fwrite(&ROW_OFFSET(matBlockDiag, 0), sizeof(BfSize), NUM_ROW_BLOCKS(matBlockDiag) + 1, fp);
  fwrite(&COL_OFFSET(matBlockDiag, 0), sizeof(BfSize), NUM_COL_BLOCKS(matBlockDiag) + 1, fp);

  /* Recursively write each diagonal blocks: */
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockDiag, k);
    bfMatDump(block, fp);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize bfMatBlockDiagGetNumRows(BfMatBlockDiag const *matBlockDiag) {
  return ROW_OFFSET(matBlockDiag, NUM_ROW_BLOCKS(matBlockDiag));
}

BfSize bfMatBlockDiagGetNumCols(BfMatBlockDiag const *matBlockDiag) {
  return COL_OFFSET(matBlockDiag, NUM_COL_BLOCKS(matBlockDiag));
}

BfMat *bfMatBlockDiagGetRowRangeCopy(BfMatBlockDiag const *matBlockDiag, BfSize i0, BfSize i1) {
  BF_ERROR_BEGIN();

  if (i0 > i1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatBlockCoo *blockRow = NULL;

  BfPtrArray indexedRowBlocks;
  bfInitPtrArrayWithDefaultCapacity(&indexedRowBlocks);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfSize i0_ = ROW_OFFSET(matBlockDiag, k);
    BfSize i1_ = ROW_OFFSET(matBlockDiag, k + 1);

    /* If the current block's row range doesn't overlap with the row
     * range to copy, skip it */
    if (i1_ <= i0 || i1 <= i0_)
      continue;

    BfSize j0_ = COL_OFFSET(matBlockDiag, k);
    BfSize m_ = NUM_BLOCK_ROWS(matBlockDiag, k);
    BfMat const *block = BLOCK(matBlockDiag, k);

    /* The range of rows which should be copied from `block` with
     * respect to `block`'s index system. */
    BfSize i0__ = i0_ < i0 ? i0 - i0_ : 0;
    BfSize i1__ = m_ - (i1 < i1_ ? i1_ - i1 : 0);

    /* Copy the row range---the whole point here is that this can
     * easily be a subset of the rows of `block`! */
    BfMat *blockRowRange = bfMatGetRowRangeCopy(block, i0__, i1__);
    HANDLE_ERROR();

    // TODO: explain how the indexing for i0 works
    BfIndexedMat *indexedRowBlock = bfIndexedMatNewFromMat(
      i0_ < i0 ? 0 : i0_ - i0,
      j0_,
      blockRowRange,
      BF_POLICY_STEAL);

    bfPtrArrayAppend(&indexedRowBlocks, indexedRowBlock);
    HANDLE_ERROR();

    bfMatDelete(&blockRowRange);
  }

  BfSize m = i1 - i0;
  BfSize n = bfMatGetNumCols(bfMatBlockDiagConstToMatConst(matBlockDiag));

  blockRow = bfMatBlockCooNewFromIndexedBlocks(m, n, &indexedRowBlocks, BF_POLICY_STEAL);
  HANDLE_ERROR();

  BF_ERROR_END() {
    bfMatBlockCooDeinitAndDealloc(&blockRow);

    BF_DIE(); // ???
  }

  for (BfSize k = 0; k < bfPtrArraySize(&indexedRowBlocks); ++k) {
    BfIndexedMat *indexedRowBlock = bfPtrArrayGet(&indexedRowBlocks, k);
    BF_ASSERT(bfMatIsView(indexedRowBlock->mat));
    bfMatDelete(&indexedRowBlock->mat);
    bfMemFree(indexedRowBlock);
  }

  bfPtrArrayDeinit(&indexedRowBlocks);

  return bfMatBlockCooToMat(blockRow);
}

void bfMatBlockDiagScaleCols(BfMat *mat, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMatBlock *matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  BfMatBlockDiag *matBlockDiag = bfMatToMatBlockDiag(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfSize j0 = matBlock->colOffset[k];
    BfSize j1 = matBlock->colOffset[k + 1];
    BfMat *block = matBlock->block[k];
    BfVec *subvec = bfVecGetSubvecCopy(vec, j0, j1);
    bfMatScaleCols(block, subvec);
    bfVecDelete(&subvec);
  }

  BF_ERROR_END() {}
}

BfMat *bfMatBlockDiagMul(BfMat const *mat, BfMat const *other) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = bfMatConstToMatBlockConst(mat);
  BfMatBlockDiag const *matBlockDiag = bfMatConstToMatBlockDiagConst(mat);

  BfSize numRows = bfMatGetNumRows(mat);
  BfSize numCols = bfMatGetNumCols(other);
  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);

  BfMat *result = NULL;
  BfMat *op2Rows = NULL;
  BfMat *resultRows = NULL;

  result = bfMatEmptyLike(other, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize i = 0, i0, i1, j0, j1; i < numBlocks; ++i) {
    i0 = matBlock->rowOffset[i];
    i1 = matBlock->rowOffset[i + 1];
    j0 = matBlock->colOffset[i];
    j1 = matBlock->colOffset[i + 1];

    op2Rows = bfMatGetRowRange((BfMat *)other, j0, j1);
    resultRows = bfMatMul(matBlock->block[i], op2Rows);
    bfMatSetRowRange(result, i0, i1, resultRows);

    bfMatDelete(&resultRows);
    bfMatDelete(&op2Rows);
  }

  BF_ERROR_END()
    bfMatDelete(&result);

  return result;
}

BfVec *bfMatBlockDiagMulVec(BfMatBlockDiag const *matBlockDiag, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatBlockDiagConstToMatConst(matBlockDiag);

  if (bfMatGetNumCols(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize resultSize = bfMatGetNumRows(mat);

  BfVec *result = NULL;
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
#if BF_DEBUG
    result = bfVecRealToVec(bfVecRealNewWithValue(resultSize, BF_NAN));
#else
    result = bfVecRealToVec(bfVecRealNewEmpty(resultSize));
#endif
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockDiag, k);

    BfSize i0 = ROW_OFFSET(matBlockDiag, k);
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = COL_OFFSET(matBlockDiag, k);
    BfSize j1 = j0 + bfMatGetNumCols(block);

    /* Multiply the current block with the relevant part of `vec`: */
    BfVec const *subvecView = bfVecGetSubvecViewConst(vec, j0, j1);
    BfVec *tmp = bfMatMulVec(block, subvecView);
    bfVecDelete((BfVec **)&subvecView);

    /* Assign the result to the relevant part of `result`: */
    bfVecSetRange(result, i0, i1, tmp);

    bfVecDelete(&tmp);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

BfVec *bfMatBlockDiagRmulVec(BfMatBlockDiag const *matBlockDiag, BfVec const *vec) {
  BF_ERROR_BEGIN();

  BfMat const *mat = bfMatBlockDiagConstToMatConst(matBlockDiag);

  if (bfMatGetNumRows(mat) != vec->size)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize resultSize = bfMatGetNumCols(mat);

  BfVec *result = NULL;
  switch (bfVecGetType(vec)) {
  case BF_TYPE_VEC_REAL:
#if BF_DEBUG
    result = bfVecRealToVec(bfVecRealNewWithValue(resultSize, BF_NAN));
#else
    result = bfVecRealToVec(bfVecRealNewEmpty(resultSize));
#endif
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockDiag, k);

    BfSize i0 = ROW_OFFSET(matBlockDiag, k);
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = COL_OFFSET(matBlockDiag, k);
    BfSize j1 = j0 + bfMatGetNumCols(block);

    /* Multiply the current block with the relevant part of `vec`: */
    BfVec const *subvecView = bfVecGetSubvecViewConst(vec, i0, i1);
    BfVec *tmp = bfMatRmulVec(block, subvecView);
    bfVecDelete((BfVec **)&subvecView);

    /* Assign the result to the relevant part of `result`: */
    bfVecSetRange(result, j0, j1, tmp);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

void bfMatBlockDiagNegate(BfMat *mat) {
  BF_ERROR_BEGIN();

  BfMatBlock *matBlock = NULL;
  BfMatBlockDiag *matBlockDiag = NULL;

  matBlock = bfMatToMatBlock(mat);
  HANDLE_ERROR();

  BfSize numBlocks = bfMatBlockNumBlocks(matBlock);

  matBlockDiag = bfMatToMatBlockDiag(mat);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numBlocks; ++i)
    bfMatNegate(bfMatBlockDiagGetBlock(matBlockDiag, i));

  BF_ERROR_END() {}
}

void bfMatBlockDiagPrintBlocksDeep(BfMatBlockDiag const *matBlockDiag, FILE *fp,
                                   BfSize i0, BfSize j0, BfSize depth) {
  BfMat const *mat = bfMatBlockDiagConstToMatConst(matBlockDiag);

  BfSize i1 = i0 + bfMatGetNumRows(mat);
  BfSize j1 = j0 + bfMatGetNumCols(mat);

  fprintf(fp, "%u %lu %lu %lu %lu %lu\n", BF_TYPE_MAT_BLOCK_DIAG, i0, i1, j0, j1, depth);

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);
  for (BfSize k = 0; k < numBlocks; ++k) {
    BfSize i0_ = i0 + ROW_OFFSET(matBlockDiag, k);
    BfSize j0_ = j0 + COL_OFFSET(matBlockDiag, k);
    BfMat const *block = matBlockDiag->super.block[k];
    bfMatPrintBlocksDeep(block, fp, i0_, j0_, depth + 1);
  }
}

BfMat *solve_matDenseComplex(BfMatBlockDiag const *matBlockDiag,
                             BfMatDenseComplex const *otherMatDenseComplex) {
  BF_ERROR_BEGIN();

  // (m x n) (n x p) = (m x p)

  BfMat const *mat = bfMatBlockDiagConstToMatConst(matBlockDiag);
  BfMat const *otherMat = bfMatDenseComplexConstToMatConst(otherMatDenseComplex);

  if (bfMatGetNumCols(mat) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize m = bfMatGetNumRows(mat);
  BfSize p = bfMatGetNumCols(otherMat);

  BfMatDenseComplex *resultMatDenseComplex = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(resultMatDenseComplex, m, p);
  HANDLE_ERROR();

  BfMat *resultMat = bfMatDenseComplexToMat(resultMatDenseComplex);

  BfSize numBlocks = bfMatBlockDiagNumBlocks(matBlockDiag);

  for (BfSize k = 0; k < numBlocks; ++k) {
    BfMat const *block = BLOCK(matBlockDiag, k);

    BfSize i0 = ROW_OFFSET(matBlockDiag, k);
    BfSize i1 = i0 + bfMatGetNumRows(block);

    BfSize j0 = COL_OFFSET(matBlockDiag, k);
    BfSize j1 = j0 + bfMatGetNumCols(block);

    BfMat const *otherBlock = bfMatGetRowRange((BfMat *)otherMat, j0, j1);
    BfMat *resultBlock = bfMatSolve(block, otherBlock);
    bfMatSetRowRange(resultMat, i0, i1, resultBlock);
    bfMatDelete(&resultBlock);
    bfMatDelete((BfMat **)&otherBlock);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return resultMat;
}

BfMat *bfMatBlockDiagSolve(BfMatBlockDiag const *matBlockDiag, BfMat const *mat) {
  switch (bfMatGetType(mat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return solve_matDenseComplex(matBlockDiag, bfMatConstToMatDenseComplexConst(mat));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

void bfMatBlockDiagTranspose(BfMatBlockDiag *matBlockDiag) {
  BF_ERROR_BEGIN();

  if (bfMatIsView(bfMatBlockDiagToMat(matBlockDiag)))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Swap the number of blocks in each row and column: */
  SWAP(matBlockDiag->super.super.numRows, matBlockDiag->super.super.numCols);

  /* Swap the row and column offsets: */
  SWAP(matBlockDiag->super.rowOffset, matBlockDiag->super.colOffset);

  /* Transpose each diagonal block: */
  for (BfSize k = 0; k < bfMatBlockDiagNumBlocks(matBlockDiag); ++k) {
    bfMatTranspose(BLOCK(matBlockDiag, k));
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

/** Interface: MatBlock */

static BfMatBlockVtable MAT_BLOCK_VTABLE = {
  .NumBlocks = (__typeof__(&bfMatBlockNumBlocks))bfMatBlockDiagNumBlocks,
  .GetNumRowBlocks = (__typeof__(&bfMatBlockGetNumRowBlocks))bfMatBlockDiagGetNumRowBlocks,
  .GetNumColBlocks = (__typeof__(&bfMatBlockGetNumColBlocks))bfMatBlockDiagGetNumColBlocks,
  .GetBlockCopy = (__typeof__(&bfMatBlockGetBlockCopy))bfMatBlockDiagGetBlockCopy,
};

BfSize bfMatBlockDiagNumBlocks(BfMatBlockDiag const *matBlockDiag) {
  BfMat const *mat = bfMatBlockDiagConstToMatConst(matBlockDiag);
  return mat->numRows < mat->numCols ? mat->numRows : mat->numCols;
}

BfSize bfMatBlockDiagGetNumRowBlocks(BfMatBlockDiag const *matBlockDiag) {
  return NUM_ROW_BLOCKS(matBlockDiag);
}

BfSize bfMatBlockDiagGetNumColBlocks(BfMatBlockDiag const *matBlockDiag) {
  return NUM_COL_BLOCKS(matBlockDiag);
}

BfMat *bfMatBlockDiagGetBlockCopy(BfMatBlockDiag const *matBlockDiag, BfSize i, BfSize j) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = bfMatBlockDiagConstToMatBlockConst(matBlockDiag);
  HANDLE_ERROR();

  BfMat *block = NULL;

  BfSize m = bfMatBlockGetNumBlockRows(matBlock, i);
  BfSize n = bfMatBlockGetNumBlockCols(matBlock, j);

  /* Copy the ith block if this is a diagonal block. */
  if (i == j) {
    if (matBlock->block[i] == NULL)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

    block = bfMatCopy(matBlock->block[i]);
    HANDLE_ERROR();
  }

  /* Otherwise, return a new zero block. */
  else {
    BfMatZero *matZero = bfMatZeroNew();
    HANDLE_ERROR();

    bfMatZeroInit(matZero, m, n);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    bfMatDelete(&block);
  }

  return block;
}

/** Upcasting: MatBlockDiag -> Mat */

BfMat *bfMatBlockDiagToMat(BfMatBlockDiag *matBlockDiag) {
  return &matBlockDiag->super.super;
}

BfMat const *bfMatBlockDiagConstToMatConst(BfMatBlockDiag const *matBlockDiag) {
  return &matBlockDiag->super.super;
}

/** Upcasting: MatBlockDiag -> MatBlock */

BfMatBlock *bfMatBlockDiagToMatBlock(BfMatBlockDiag *matBlockDiag) {
  return &matBlockDiag->super;
}

BfMatBlock const *bfMatBlockDiagConstToMatBlockConst(BfMatBlockDiag const *matBlockDiag) {
  return &matBlockDiag->super;
}

/** Downcasting: Mat -> MatBlockDiag */

BfMatBlockDiag *bfMatToMatBlockDiag(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_DIAG)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDiag *)mat;
  }
}

BfMatBlockDiag const *bfMatConstToMatBlockDiagConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_BLOCK_DIAG)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatBlockDiag const *)mat;
  }
}

/** Implementation: MatBlockDiag */

BfMatBlockDiag *bfMatBlockDiagNew() {
  BF_ERROR_BEGIN();

  BfMatBlockDiag *mat = bfMemAlloc(1, sizeof(BfMatBlockDiag));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return mat;
}

BfMatBlockDiag *bfMatBlockDiagNewFromBlocks(BfPtrArray *blocks, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfMatBlockDiag *matBlockDiag = bfMatBlockDiagNew();
  HANDLE_ERROR();

  BfSize numBlocks = bfPtrArraySize(blocks);

  for (BfSize i = 0; i < numBlocks; ++i) {
    BfMat const *block = bfPtrArrayGet(blocks, i);
    if (block == NULL)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  bfMatBlockInit(&matBlockDiag->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE,
                 numBlocks, numBlocks, numBlocks);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numBlocks; ++i) {
    BfMat *block = bfPtrArrayGet(blocks, i);
    BLOCK(matBlockDiag, i) = bfMatGet(block, policy);
  }

  ROW_OFFSET(matBlockDiag, 0) = 0;
  COL_OFFSET(matBlockDiag, 0) = 0;
  for (BfSize i = 0; i < numBlocks; ++i) {
    BfMat const *block = bfPtrArrayGet(blocks, i);
    ROW_OFFSET(matBlockDiag, i + 1) = bfMatGetNumRows(block);
    COL_OFFSET(matBlockDiag, i + 1) = bfMatGetNumCols(block);
  }
  bfSizeRunningSum(numBlocks + 1, &ROW_OFFSET(matBlockDiag, 0));
  bfSizeRunningSum(numBlocks + 1, &COL_OFFSET(matBlockDiag, 0));

  BF_ERROR_END() {
    BF_DIE();
  }

  return matBlockDiag;
}

void bfMatBlockDiagInit(BfMatBlockDiag *mat, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfSize numBlocks = numRows < numCols ? numRows : numCols;

  bfMatBlockInit(&mat->super, &MAT_VTABLE, &MAT_BLOCK_VTABLE, numBlocks, numRows, numCols);
  HANDLE_ERROR();

  BF_ERROR_END() {}
}

void bfMatBlockDiagDeinit(BfMatBlockDiag *mat) {
  bfMatBlockDeinit(&mat->super);
}

void bfMatBlockDiagDealloc(BfMatBlockDiag **mat) {
  bfMemFree(*mat);
  *mat = NULL;
}

void bfMatBlockDiagDeinitAndDealloc(BfMatBlockDiag **mat) {
  bfMatBlockDiagDeinit(*mat);
  bfMatBlockDiagDealloc(mat);
}

BfMat *bfMatBlockDiagGetBlock(BfMatBlockDiag *matBlockDiag, BfSize i) {
  BF_ERROR_BEGIN();

  BfMatBlock *matBlock = bfMatBlockDiagToMatBlock(matBlockDiag);
  BfMat *block = NULL;

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (i >= numRowBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  if (i >= numColBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  block = matBlockDiag->super.block[i];

  BF_ERROR_END()
    block = NULL;

  return block;
}

BfMat const *bfMatBlockDiagGetBlockConst(BfMatBlockDiag const *matBlockDiag, BfSize i) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = bfMatBlockDiagConstToMatBlockConst(matBlockDiag);
  BfMat const *block = NULL;

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (i >= numRowBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  if (i >= numColBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  block = bfMatGet(matBlockDiag->super.block[i], BF_POLICY_VIEW);

  BF_ERROR_END() {
    BF_DIE();
  }

  return block;
}

void bfMatBlockDiagSetBlock(BfMatBlockDiag *matBlockDiag, BfSize i, BfMat *mat) {
  BF_ERROR_BEGIN();

  BfMatBlock const *matBlock = bfMatBlockDiagConstToMatBlockConst(matBlockDiag);

  BfSize numRowBlocks = bfMatBlockGetNumRowBlocks(matBlock);
  if (i >= numRowBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numColBlocks = bfMatBlockGetNumColBlocks(matBlock);
  if (i >= numColBlocks)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  matBlockDiag->super.block[i] = mat;

  BF_ERROR_END() {}
}
