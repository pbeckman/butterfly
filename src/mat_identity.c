#include <bf/mat_identity.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .Delete = (__typeof__(&bfMatDelete))bfMatIdentityDelete,
  .GetType = (__typeof__(&bfMatGetType))bfMatIdentityGetType,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatIdentityGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatIdentityGetNumCols,
  .GetRowRangeCopy = (__typeof__(&bfMatGetRowRangeCopy))bfMatIdentityGetRowRangeCopy,
};

void bfMatIdentityDelete(BfMatIdentity **matIdentity) {
  bfMatIdentityDeinitAndDealloc(matIdentity);
}

BfType bfMatIdentityGetType(BfMatIdentity const *matIdentity) {
  (void)matIdentity;
  return BF_TYPE_MAT_IDENTITY;
}

BfSize bfMatIdentityGetNumRows(BfMatIdentity const *matIdentity) {
  return matIdentity->super.numRows;
}

BfSize bfMatIdentityGetNumCols(BfMatIdentity const *matIdentity) {
  return matIdentity->super.numCols;
}

BfMat *bfMatIdentityGetRowRangeCopy(BfMatIdentity const *matIdentity, BfSize i0, BfSize i1) {
  BEGIN_ERROR_HANDLING();

  BfMat const *mat = bfMatIdentityConstToMatConst(matIdentity);

  BfMat *rowRangeCopy = NULL;

  BfSize m = bfMatGetNumRows(mat);

  if (i0 == 0 && i1 == m) {
    BfMatIdentity *matIdentity_ = bfMatIdentityNew();
    HANDLE_ERROR();

    bfMatIdentityInit(matIdentity_, m);
    HANDLE_ERROR();

    rowRangeCopy = bfMatIdentityToMat(matIdentity_);
  } else {
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}

  return rowRangeCopy;
}

/** Upcasting: MatIdentity -> Mat */

BfMat *bfMatIdentityToMat(BfMatIdentity *matIdentity) {
  return &matIdentity->super;
}

BfMat const *bfMatIdentityConstToMatConst(BfMatIdentity const *matIdentity) {
  return &matIdentity->super;
}

/** Downcasting: Mat -> MatIdentity */

BfMatIdentity *bfMatToMatIdentity(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_IDENTITY)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfMatIdentity *)mat;
  }
}

/** Implementation: MatIdentity */

BfMatIdentity *bfMatIdentityNew() {
  BEGIN_ERROR_HANDLING();

  BfMatIdentity *matIdentity = malloc(sizeof(BfMatIdentity));
  if (matIdentity == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return matIdentity;
}

void bfMatIdentityInit(BfMatIdentity *matIdentity, BfSize n) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&matIdentity->super, &MAT_VTABLE, n, n);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfMatIdentityDeinit(matIdentity);
  }
}

void bfMatIdentityDeinit(BfMatIdentity *matIdentity) {
  (void)matIdentity;
  assert(false);
}

void bfMatIdentityDealloc(BfMatIdentity **matIdentity) {
  free(*matIdentity);
  *matIdentity = NULL;
}

void bfMatIdentityDeinitAndDealloc(BfMatIdentity **matIdentity) {
  bfMatIdentityDeinit(*matIdentity);
  bfMatIdentityDealloc(matIdentity);
}
