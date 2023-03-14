#pragma once

#include "mat.h"
#include "tree.h"

typedef struct BfFacStreamer BfFacStreamer;

typedef struct BfFacStreamerSpec {
  BfTree *rowTree;
  BfTree *colTree;
  BfSize rowTreeInitDepth;
  BfSize colTreeInitDepth;
  BfReal tol;
  BfSize minNumRows;
  BfSize minNumCols;
} BfFacStreamerSpec;

BfFacStreamer *bfFacStreamerNew();
void bfFacStreamerInit(BfFacStreamer *facStreamer, BfFacStreamerSpec const *spec);
void bfFacStreamerDeinit(BfFacStreamer *facStreamer);
BfSize bfFacStreamerGetNumRows(BfFacStreamer const *facStreamer);
void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat const *Phi);
bool bfFacStreamerIsDone(BfFacStreamer const *facStreamer);
BfMat *bfFacStreamerGetFac(BfFacStreamer const *facStreamer);
BfTreeNode *bfFacStreamerGetCurrentColumnNode(BfFacStreamer const *facStreamer);
BfPerm const *bfFacStreamerGetRowTreeReversePerm(BfFacStreamer const *facStreamer);
