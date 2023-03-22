#include <bf/eval_tree.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <bf/cheb.h>
#include <bf/error.h>
#include <bf/error_macros.h>

// NOTE: a better way to do these is to allocate a pool of coefficient
// vectors and BfEvalTreeNodes in BfEvalTree... can use bitfields to
// combine "isLeaf" and the index into one. The idea here is that
// since every internal node has the same number of children, we only
// need the index to the starting position.

typedef struct BfEvalTreeNode BfEvalTreeNode;

struct BfEvalTreeNode {
  BfReal a, b;
  bool isLeaf;
  BfChebStd *cheb;
  BfEvalTreeNode *children;
};

struct BfEvalTree {
  BfReal (*f)(BfReal);
  BfReal a, b;
  BfSize d;
  BfSize k;
  BfReal tol;
  BfReal *x; // chebyshev points

  BfEvalTreeNode *root;
};

void evalTreeNodeInitRec(BfEvalTreeNode *node, BfEvalTree const *tree) {
  node->cheb = bfChebStdNewWithDegree(tree->d);
  bfChebStdInterp(node->cheb, tree->f, node->a, node->b, tree->x);

  /* Check the interpolation error. If it's small enough, declare this
   * a leaf node and return. */
  BfReal interpError = bfChebStdGetErrorEstimate(node->cheb);
  node->isLeaf = interpError <= tree->tol;
  if (node->isLeaf) {
    node->children = NULL;
    return;
  }

  /** Otherwise, we need to free the current Chebyshev polynomial and
   * recursively create the cchildren. */

  bfChebStdDelete(&node->cheb);

  node->children = malloc(tree->k*sizeof(BfEvalTreeNode));

  BfReal delta = node->b - node->a;
  delta /= tree->k;
  for (BfSize i = 0; i < tree->k; ++i) {
    node->children[i].a = i == 0 ? node->a : node->a + i*delta;
    node->children[i].b = i == tree->k - 1 ? node->b : node->a + (i + 1)*delta;
    evalTreeNodeInitRec(&node->children[i], tree);
  }

  if (node->isLeaf) {
    assert(node->cheb != NULL);
    assert(node->children == NULL);
  } else {
    assert(node->cheb == NULL);
    assert(node->children != NULL);
  }
}

BfReal evalTreeNodeGetValueRec(BfEvalTreeNode const *node, BfEvalTree const *tree, BfReal x) {
  if (node->isLeaf) {
    assert(node->a <= x && x <= node->b);
    return bfChebStdEval(node->cheb, 2*(x - node->a)/(node->b - node->a) - 1);
  } else {
    for (BfSize i = 0; i < tree->k; ++i) {
      BfEvalTreeNode const *child = &node->children[i];
      if (child->a <= x && x <= child->b)
        return evalTreeNodeGetValueRec(child, tree, x);
    }
  }
  assert(false);
}

BfReal bfEvalTreeGetValue(BfEvalTree const *tree, BfReal x) {
  return evalTreeNodeGetValueRec(tree->root, tree, x);
}

BfEvalTree *bfEvalTreeNew() {
  BEGIN_ERROR_HANDLING();

  BfEvalTree *evalTree = malloc(sizeof(BfEvalTree));
  if (evalTree == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {
    assert(false);
  }

  return evalTree;
}

void bfEvalTreeInit(BfEvalTree *tree, BfEvalTreeSpec const *spec) {
  tree->f = spec->f;
  tree->a = spec->a;
  tree->b = spec->b;
  tree->d = spec->d;
  tree->k = spec->k;
  tree->tol = spec->tol;

  tree->x = malloc((tree->d + 1)*sizeof(BfReal));
  bfGetChebPts(tree->d + 1, tree->x);

  tree->root = malloc(sizeof(BfEvalTreeNode));
  tree->root->a = tree->a;
  tree->root->b = tree->b;
  evalTreeNodeInitRec(tree->root, tree);
}