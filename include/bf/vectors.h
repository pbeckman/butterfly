#pragma once

#include "geom.h"
#include "mat_dense_real.h"

void bfVector2Normalize(BfVector2 u);
void bfVector2Reject(BfVector2 u, BfVector2 const v);
void bfVector2Negate(BfVector2 u);
void bfVector2Rotate(BfVector2 u, BfReal theta);

void bfVector3Copy(BfVector3 u, BfVector3 const v);
void bfVector3Scale(BfVector3 u, BfReal alpha);
BfReal bfVector3Norm(BfVector3 const u);
BfReal bfVector3Dot(BfVector3 const u, BfVector3 const v);
void bfVector3Cross(BfVector3 const u, BfVector3 const v, BfVector3 w);
void bfVector3Normalize(BfVector3 u);
void bfVector3Negate(BfVector3 u);

struct BfVectors2 {
  BfVector2 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfVectors2 *bfVectors2NewEmpty(void);
BfVectors2 *bfVectors2NewWithCapacity(BfSize capacity);
BfVectors2 *bfVectors2NewFromFile(char const *path);
BfVectors2 const *bfVectors2ConstViewFromMat(BfMat const *mat);
BfVectors2 const *bfVectors2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal);
BfVectors2 bfGetUninitializedVectors2(void);
void bfVectors2InitEmpty(BfVectors2 *vectors);
void bfVectors2InitWithCapacity(BfVectors2 *vectors, BfSize capacity);
void bfVectors2InitFromFile(char const *path, BfVectors2 *vectors);
void bfVectors2Deinit(BfVectors2 *vectors);
void bfVectors2Dealloc(BfVectors2 **vectors);
void bfVectors2DeinitAndDealloc(BfVectors2 **vectors);
void bfGetVectorsByIndex(BfVectors2 const *vectors,
                        BfSize numInds, BfSize const *inds,
                        BfVectors2 *indexedVectors);
void bfSaveVectors2(BfVectors2 const *vectors, char const *path);
void bfVectors2Append(BfVectors2 *vectors, BfVector2 const v);
void bfVectors2Extend(BfVectors2 *vectors, BfVectors2 const *newVectors);
void bfVectors2Get(BfVectors2 const *vectors, BfSize i, BfVector2 v);
BfSize bfVectors2GetSize(BfVectors2 const *vectors);
void bfVectors2Set(BfVectors2 *vectors, BfSize i, BfVector2 const v);
BfVectors2 *bfVectors2GetRangeView(BfVectors2 *vectors, BfSize i0, BfSize i1);

struct BfVectors3 {
  BfVector3 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfVectors3 *bfVectors3NewEmpty(void);
BfVectors3 *bfVectors3NewWithCapacity(BfSize capacity);
BfVectors3 *bfVectors3Copy(BfVectors3 const *vectors);
void bfVectors3InitEmpty(BfVectors3 *vectors);
void bfVectors3InitWithCapacity(BfVectors3 *vectors, BfSize capacity);
void bfVectors3Deinit(BfVectors3 *vectors);
void bfVectors3Dealloc(BfVectors3 **vectors);
void bfVectors3DeinitAndDealloc(BfVectors3 **vectors);
void bfVectors3GetByIndex(BfVectors3 const *vectors, BfSize numInds, BfSize const *inds, BfVectors3 *indexedVectors);
void bfVectors3Append(BfVectors3 *vectors, BfVector3 const v);
BfReal const *bfVectors3GetConstPtr(BfVectors3 const *vectors, BfSize i);
