#include <bf/assert.h>
#include <bf/error.h>
#include <bf/fiedler_tree.h>
#include <bf/fiedler_tree_node.h>
#include <bf/logging.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/util.h>

#include <stdlib.h>

#define DEBUG_OUTPUT 0

void checkPerm(BfTree const *tree, BfTreeNode const *treeNode, void *arg);
void writeNodeSpan(BfTree const *tree, BfTreeNode const *treeNode, void *arg);
void incrementNumLevelNodes(BfTree const *tree, BfTreeNode const *treeNode, void *arg);

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 3) {
    printf("usage: %s <trimesh.obj> [tol]\n", argv[0]);
    return 1;
  }

  bfSetLogLevel(BF_LOG_LEVEL_DEBUG);

  char const *objPath = argv[1];

  BfReal tol = argc >= 3 ? atof(argv[2]) : 1e-15;

  printf("running Fiedler tree test:\n");

  bfToc();

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(objPath);
  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  printf("- loaded triangle mesh (%lu verts and %lu faces) [%.1fs]\n",
         numVerts, bfTrimeshGetNumFaces(trimesh), bfToc());

  BfFiedlerTree *fiedlerTree = bfFiedlerTreeNewFromTrimesh(trimesh, tol, /* keepNodeTrimeshes: */ true);
  printf("- built Fiedler tree (tol = %g) [%.1fs]\n", tol, bfToc());

  BfTree *tree = bfFiedlerTreeToTree(fiedlerTree);
  bfTreeMapConst(tree, NULL, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfTreeMapConstFunc)checkPerm, &tol);

  BfSize maxDepth = bfTreeGetMaxDepth(tree);
  printf("- depth of deepest node: %lu\n", maxDepth);

  BfSize *numLevelNodes = bfMemAllocAndZero(maxDepth + 1, sizeof(BfSize));

  bfTreeMapConst(tree, NULL, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfTreeMapConstFunc)incrementNumLevelNodes, numLevelNodes);

  BF_ASSERT(numLevelNodes[0] == numVerts);
  for (BfSize depth = 0; depth < maxDepth; ++depth)
    BF_ASSERT(numLevelNodes[depth + 1] <= numLevelNodes[depth]);

  BfSize maxDepthComplete = 0;
  for (BfSize depth = 1; depth <= maxDepth; ++depth)
    if (numLevelNodes[depth] == numVerts)
      ++maxDepthComplete;

  printf("- deepest complete level: %lu\n", maxDepthComplete);

  bfMemFree(numLevelNodes);

  FILE *fp = fopen("node_spans.txt", "w");
  bfTreeMapConst(tree, NULL, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfTreeMapConstFunc)writeNodeSpan, fp);
  fclose(fp);

  bfPermSave(bfTreeGetPermConst(tree), "perm.bin");

  bfFiedlerTreeDeinitAndDealloc(&fiedlerTree);
  bfTrimeshDeinitAndDealloc(&trimesh);
}

void checkPerm(BfTree const *tree, BfTreeNode const *treeNode, void *arg) {
  BfReal tol = *(BfReal *)arg;

  BfFiedlerTree const *fiedlerTree = bfTreeConstToFiedlerTreeConst(tree);
  BfFiedlerTreeNode const *fiedlerTreeNode = bfTreeNodeConstToFiedlerTreeNodeConst(treeNode);

  BfSize i0 = bfTreeNodeGetFirstIndex(treeNode);
  BfSize i1 = bfTreeNodeGetLastIndex(treeNode);

  for (BfSize i = i0; i < i1; ++i) {
    BfReal const *v = bfTrimeshGetVertPtrConst(fiedlerTree->trimesh, tree->perm->index[i]);
    BF_ASSERT(bfPoints3ContainsApprox(bfTrimeshGetVertsConst(fiedlerTreeNode->trimesh), v, tol));
  }
}

void writeNodeSpan(BfTree const *tree, BfTreeNode const *treeNode, void *arg) {
  BfSize depth = bfTreeNodeGetDepth(treeNode);
  BfSize i0 = bfTreeNodeGetFirstIndex(treeNode);
  BfSize i1 = bfTreeNodeGetLastIndex(treeNode);
  fprintf((FILE *)arg, "%lu %lu %lu\n", depth, i0, i1);
}

void incrementNumLevelNodes(BfTree const *tree, BfTreeNode const *treeNode, void *arg) {
  BfSize *numLevelNodes = (BfSize *)arg;
  BfSize depth = bfTreeNodeGetDepth(treeNode);
  numLevelNodes[depth] += bfTreeNodeGetNumPoints(treeNode);
}
