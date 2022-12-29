#include <bf/interval_tree.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/interval_tree_node.h>
#include <bf/tree.h>

/** Interface: Tree -> IntervalTree */

static BfTreeVtable TreeVtable = {
  .GetType = (__typeof__(&bfTreeGetType))bfIntervalTreeGetType
};

BfType bfIntervalTreeGetType(BfTree const *tree) {
  (void)tree;
  return BF_TYPE_INTERVAL_TREE;
}

/** Upcasting: */

BfTree *bfIntervalTreeToTree(BfIntervalTree *intervalTree) {
  return &intervalTree->super;
}

/** Implementation: IntervalTree */

BfIntervalTree *bfIntervalTreeNew() {
  BEGIN_ERROR_HANDLING();

  BfIntervalTree *intervalTree = malloc(sizeof(BfIntervalTree));
  if (intervalTree == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    intervalTree = NULL;

  return intervalTree;
}

void bfIntervalTreeInitEmpty(BfIntervalTree *intervalTree, BfReal a, BfReal b, BfSize k) {
  BEGIN_ERROR_HANDLING();

  BfIntervalTreeNode *root = bfIntervalTreeNodeNew();
  HANDLE_ERROR();

  bfIntervalTreeNodeInitEmpty(root, intervalTree, a, b, k);

  bfTreeInit(&intervalTree->super, &TreeVtable, bfIntervalTreeNodeToTreeNode(root), 0);

  END_ERROR_HANDLING() {}
}
