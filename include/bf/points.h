#pragma once

#include "geom.h"
#include "mat.h"

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);
BfReal bfPoint2Magnitude(BfPoint2 const p);
void bfPoint2SampleUniformlyFromBoundingBox(BfBbox2 const *bbox, BfPoint2 p);
void bfPoint2RotateAboutOrigin(BfPoint2 p, BfReal theta);
void bfPoint2Translate(BfPoint2 p, BfVector2 const u);

BfReal bfPoint3Dist(BfPoint3 const p, BfPoint3 const q);
void bfPoint3Sub(BfPoint3 const v, BfPoint3 const u, BfVector3 uv);
void bfPoint3GetPointOnRay(BfPoint3 const r0, BfVector3 const dr, BfReal t, BfPoint3 rt);
void bfPoint3Copy(BfPoint3 x, BfPoint3 const y);

struct BfPoints1 {
  BfPoint1 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfPoints1 *bfPoints1New();
void bfPoints1InitEmpty(BfPoints1 *points, BfSize capacity);
void bfPoints1InitViewFromVecReal(BfPoints1 *points, BfVecReal const *vecReal);
void bfPoints1Deinit(BfPoints1 *points);
bool bfPoints1IsSorted(BfPoints1 const *points);
void bfPoints1Append(BfPoints1 *points, BfPoint1 point);
void bfPoints1InsertPointsSorted(BfPoints1 *points, BfPoints1 const *newPoints);

struct BfPoints2 {
  BfPoint2 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfPoints2 *bfPoints2NewEmpty();
BfPoints2 *bfPoints2NewGrid(BfBbox2 const *bbox, BfSize nx, BfSize ny);
BfPoints2 const *bfPoints2ConstViewFromMat(BfMat const *mat);
BfPoints2 const *bfPoints2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal);
BfPoints2 bfGetUninitializedPoints2();
void bfInitEmptyPoints2(BfPoints2 *points, BfSize numPoints);
void bfReadPoints2FromFile(char const *path, BfPoints2 *points);
void bfFreePoints2(BfPoints2 *points);
bool bfPoints2Initialized(BfPoints2 const *points);
BfBbox2 bfPoints2GetBoundingBox(BfPoints2 const *points);
void bfGetPointsByIndex(BfPoints2 const *points, BfSize numInds, BfSize const *inds, BfPoints2 *indexedPoints);
void bfPrintPoints2(BfPoints2 const *points);
void bfSavePoints2(BfPoints2 const *points, char const *path);
BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y);
void bfPoints2Append(BfPoints2 *points, BfPoint2 const p);
void bfPoints2Extend(BfPoints2 *points, BfPoints2 const *newPoints);
void bfPoints2Get(BfPoints2 const *points, BfSize i, BfPoint2 p);
BfSize bfPoints2GetSize(BfPoints2 const *points);

struct BfPoints3 {
  BfPoint3 *data;
  BfSize size;
};

void bfPoints3InitEmpty(BfPoints3 *points, BfSize numPoints);
void bfPoints3InitFromBinaryFile(BfPoints3 *points, char const *path);
void bfPoints3Deinit(BfPoints3 *points);
BfBoundingBox3 bfPoints3GetBoundingBox(BfPoints3 const *points);
void bfPoints3GetByIndex(BfPoints3 const *points, BfSize numInds, BfSize const *inds, BfPoints3 *indexedPoints);
