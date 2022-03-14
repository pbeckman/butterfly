#pragma once

#include "def.h"

// TODO: rename: BF_MAT_TYPE_* -> BF_TYPE_MAT_* (so that first type is
// BF_TYPE_MAT, which is consistent with the rest of the typenames)
typedef enum BfMatTypes {
  BF_MAT_TYPE_MAT,
  BF_MAT_TYPE_DENSE_COMPLEX,
  BF_MAT_TYPE_DIAG_REAL,
  BF_MAT_TYPE_BLOCK,
  BF_MAT_TYPE_BLOCK_COO,
  BF_MAT_TYPE_BLOCK_DENSE,
  BF_MAT_TYPE_BLOCK_DIAG,
  BF_MAT_TYPE_COUNT
} BfMatType;

// TODO: this is a little screwed up. The convention should probably be:
//
//   BfMat<Structure><Type>
//
// with Type = Block (or Blk... so: "BfMatDenseBlk") if we just keep a
// pointer to BfMat
typedef struct BfMat BfMat;
typedef struct BfMatBlock BfMatBlock;
typedef struct BfMatBlockCoo BfMatBlockCoo;
typedef struct BfMatBlockDense BfMatBlockDense;
typedef struct BfMatBlockDiag BfMatBlockDiag;
typedef struct BfMatDenseComplex BfMatDenseComplex;
typedef struct BfMatDiagReal BfMatDiagReal;

typedef struct BfMatVtable {
  BfMat *(*zerosLike)(BfMat const *, BfSize, BfSize);
  void (*deinit)(BfMat *);
  void (*delete)(BfMat **);
  void (*deinitAndDelete)(BfMat **);
  BfMatType (*getType)(BfMat const *);
  BfSize (*numBytes)(BfMat const *);
  void (*save)(BfMat const *, char const *);
  BfSize (*getNumRows)(BfMat const *);
  BfSize (*getNumCols)(BfMat const *);
  BfMat *(*getRowRange)(BfMat *, BfSize, BfSize);
  BfMat *(*getColRange)(BfMat *, BfSize, BfSize);
  void (*addInplace)(BfMat *, BfMat const *);
  BfMat *(*mul)(BfMat const *, BfMat const *);
  BfMat *(*lstSq)(BfMat const *, BfMat const *);
} BfMatVtable;

#define BF_DEFINE_MAT_VTABLE(Subtype)                               \
  static BfMatVtable matVtbl = {                                        \
    .zerosLike = (__typeof__(&bfMatZerosLike))bf##Subtype##ZerosLike,   \
    .deinit = (__typeof__(&bfMatDeinit))bf##Subtype##Deinit,            \
    .delete = (__typeof__(&bfMatDelete))bf##Subtype##Delete,            \
    .deinitAndDelete = (__typeof__(&bfMatDeinitAndDelete))bf##Subtype##DeinitAndDelete, \
    .getType = (__typeof__(&bfMatGetType))bf##Subtype##GetType,         \
    .numBytes = (__typeof__(&bfMatNumBytes))bf##Subtype##NumBytes,      \
    .save = (__typeof__(&bfMatSave))bf##Subtype##Save,                  \
    .getNumRows = (__typeof__(&bfMatGetNumRows))bf##Subtype##GetNumRows, \
    .getNumCols = (__typeof__(&bfMatGetNumCols))bf##Subtype##GetNumCols, \
    .getRowRange = (__typeof__(&bfMatGetRowRange))bf##Subtype##GetRowRange, \
    .getColRange = (__typeof__(&bfMatGetColRange))bf##Subtype##GetColRange, \
    .addInplace = (__typeof__(&bfMatAddInplace))bf##Subtype##AddInplace, \
    .mul = (__typeof__(&bfMatMul))bf##Subtype##Mul,                     \
    .lstSq = (__typeof__(&bfMatLstSq))bf##Subtype##LstSq                \
  };

typedef enum BfMatProps {
  BF_MAT_PROPS_NONE = 0,
  BF_MAT_PROPS_VIEW = 1 << 0,
  BF_MAT_PROPS_TRANS = 1 << 1,
  BF_MAT_PROPS_CONJ = 1 << 2,
  BF_MAT_PROPS_ORTHO = 1 << 3
} BfMatProps;

struct BfMat {
  BfMatVtable *vtbl;
  BfMatProps props;
  BfSize numRows;
  BfSize numCols;
#if BF_DEBUG
  void *aux; /* pointer to extra user-defined data for purposes of debugging */
#endif
};

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols);
bool bfMatIsTransposed(BfMat const *mat);
BfMat *bfMatConjTrans(BfMat *mat);

/* BfMat interface: */
BfMat *bfMatZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols);
void bfMatDeinit(BfMat *mat);
void bfMatDelete(BfMat **mat);
void bfMatDeinitAndDelete(BfMat **mat);
BfMatType bfMatGetType(BfMat const *mat);
BfSize bfMatNumBytes(BfMat const *mat);
void bfMatSave(BfMat const *mat, char const *path);
BfSize bfMatGetNumRows(BfMat const *mat);
BfSize bfMatGetNumCols(BfMat const *mat);
BfMat *bfMatGetRowRange(BfMat *mat, BfSize i0, BfSize i1);
BfMat *bfMatGetColRange(BfMat *mat, BfSize j0, BfSize j1);
void bfMatAddInplace(BfMat *lhs, BfMat const *rhs);
BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs);
BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs);
