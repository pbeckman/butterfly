#pragma once

#include "def.h"

typedef struct BfEvalTree BfEvalTree;

typedef struct BfEvalTreeSpec {
  /* Function to interpolate: */
  BfReal (*f)(BfReal);

  /* Left endpoint of interpolation interval: */
  BfReal a;

  /* Right endpoint of interpolation interval: */
  BfReal b;

  /* Chebysev polynomial degree: */
  BfSize d;

  /* Valence of evaluation tree: */
  BfSize k;

  /* Tolerance used for Chebyshev approximation: */
  BfReal tol;
} BfEvalTreeSpec;

BfEvalTree *bfEvalTreeNew();
void bfEvalTreeInit(BfEvalTree *evalTree, BfEvalTreeSpec const *spec);
BfReal bfEvalTreeGetValue(BfEvalTree const *evalTree, BfReal x);