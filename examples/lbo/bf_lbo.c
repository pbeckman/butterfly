#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fac_span.h>
#include <bf/fac_streamer.h>
#include <bf/fiedler_tree.h>
#include <bf/interval_tree.h>
#include <bf/interval_tree_node.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/mat_csr_real.h>
#include <bf/octree.h>
#include <bf/rand.h>
#include <bf/util.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "argtable3.h"

int const MAX_NUM_ARG_ERRORS = 20;

typedef struct {
  BfSize seed;
  BfLogLevel logLevel;

  char const *objPath;
  char const *matPath;

  bool useOctree;
  bool useFiedlerTree;

  BfReal tol;
  BfSize numLeafNodes;
  BfSize rowTreeOffset;
  BfSize freqTreeDepth;
  BfSize freqTreeOffset;

  bool compareRelativeErrors;
  bool bailIfBottomedOut;
} Opts;

bool parseArgs(Opts *opts, int argc, char *argv[]) {
  bool success = true;

  struct arg_lit *help;

  struct arg_int *seed;
  struct arg_str *logLevel;

  struct arg_str *objPath;
  struct arg_str *matPath;

  struct arg_lit *useOctree;
  struct arg_lit *useFiedlerTree;

  struct arg_dbl *tol;
  struct arg_int *numLeafNodes;
  struct arg_int *rowTreeOffset;
  struct arg_int *freqTreeDepth;
  struct arg_int *freqTreeOffset;

  struct arg_lit *compareRelativeErrors;
  struct arg_lit *bailIfBottomedOut;

  struct arg_end *end;

  void *argtable[] = {
    help = arg_lit0(NULL, "help", "Display help and exit"),

    seed = arg_int0(NULL, "seed", NULL, "Seed for random number generator (default: 0)"),
    logLevel = arg_str0(NULL, "logLevel", NULL, NULL),

    objPath = arg_str0(NULL, "objPath", NULL, NULL),
    matPath = arg_str0(NULL, "matPath", NULL, NULL),

    useOctree = arg_lit0(NULL, "useOctree", NULL),
    useFiedlerTree = arg_lit0(NULL, "useFiedlerTree", NULL),

    tol = arg_dbl0(NULL, "tol", NULL, NULL),
    numLeafNodes = arg_int0(NULL, "numLeafNodes", NULL, NULL),
    rowTreeOffset = arg_int0(NULL, "rowTreeOffset", NULL, NULL),
    freqTreeDepth = arg_int0(NULL, "freqTreeDepth", NULL, NULL),
    freqTreeOffset = arg_int0(NULL, "freqTreeOffset", NULL, NULL),

    compareRelativeErrors = arg_lit0(NULL, "compareRelativeErrors", NULL),
    bailIfBottomedOut = arg_lit0(NULL, "bailIfBottomedOut", NULL),

    end = arg_end(MAX_NUM_ARG_ERRORS)
  };

  *seed->ival = 0;
  *logLevel->sval = "error";
  *objPath->sval = "";
  *matPath->sval = "";

  useOctree->count = 1;
  useFiedlerTree->count = 0;

  *tol->dval = 1e-3;
  *numLeafNodes->ival = -1;
  *rowTreeOffset->ival = 0;
  *freqTreeDepth->ival = -1;
  *freqTreeOffset->ival = -1;

  compareRelativeErrors->count = 0;
  bailIfBottomedOut->count = 0;

  BfSize numErrors = arg_parse(argc, argv, argtable);

  if (help->count > 0) {
    printf("Usage: %s\n", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Test driver for butterfly compression of LBO\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    success = false;
    goto cleanup;
  }

  if (numErrors > 0) {
    arg_print_errors(stdout, end, argv[0]);
    printf("Try '%s --help' for more information.\n", argv[0]);
    success = false;
    goto cleanup;
  }

  if (*freqTreeDepth->ival >= 0 && *freqTreeOffset->ival >= 0) {
    printf("Only one of --freqTreeDepth or --freqTreeOffset should be set\n");
    success = false;
    goto cleanup;
  }

  if (!(useOctree->count ^ useFiedlerTree->count)) {
    printf("Pass exactly one of --useOctree or --useFiedlerTree\n");
    success = false;
    goto cleanup;
  }

  if (!(strcmp(*objPath->sval, "") ^ strcmp(*matPath->sval, ""))) {
    printf("Exactly one of --objPath or --matPath should be set\n");
    success = false;
    goto cleanup;
  }

  if (strcmp(*matPath->sval, "") && useFiedlerTree->count) {
    printf("Argument --useFiedlerTree is only implemented with --objPath, not --matPath\n");
    success = false;
    goto cleanup;
  }

  /* If neither --freqTreeDepth nor --freqTreeOffset are set, then we
   * use a value of freqTreeOffset equal to 2 by default: */
  if (*freqTreeDepth->ival == -1 && *freqTreeOffset->ival == -1) {
    *freqTreeOffset->ival = 2;
  }

  opts->seed = *seed->ival;

  if (!strcmp(*logLevel->sval, "todo")) {
    opts->logLevel = BF_LOG_LEVEL_TODO;
  } else if (!strcmp(*logLevel->sval, "debug")) {
    opts->logLevel = BF_LOG_LEVEL_DEBUG;
  } else if (!strcmp(*logLevel->sval, "info")) {
    opts->logLevel = BF_LOG_LEVEL_INFO;
  } else if (!strcmp(*logLevel->sval, "warn")) {
    opts->logLevel = BF_LOG_LEVEL_WARN;
  } else if (!strcmp(*logLevel->sval, "error")) {
    opts->logLevel = BF_LOG_LEVEL_ERROR;
  } else {
    printf("--logLevel must be one of: \"todo\", \"debug\", \"info\", \"warn\", \"error\"\n");
    success = false;
    goto cleanup;
  }

  opts->objPath = *objPath->sval;
  opts->matPath  = *matPath->sval;

  opts->useOctree = useOctree->count > 0;
  opts->useFiedlerTree = useFiedlerTree->count > 0;

  opts->tol = *tol->dval;
  opts->numLeafNodes = *numLeafNodes->ival;
  opts->rowTreeOffset = *rowTreeOffset->ival;

  if (*freqTreeDepth->ival == -1) {
    opts->freqTreeDepth = BF_SIZE_BAD_VALUE;
  } else {
    BF_ASSERT(*freqTreeDepth->ival >= 0);
    opts->freqTreeDepth = *freqTreeDepth->ival;
  }

  if (*freqTreeOffset->ival == -1) {
    opts->freqTreeOffset = BF_SIZE_BAD_VALUE;
  } else {
    BF_ASSERT(*freqTreeOffset->ival >= 0);
    opts->freqTreeOffset = *freqTreeOffset->ival;
  }

  opts->compareRelativeErrors = compareRelativeErrors->count > 0;
  opts->bailIfBottomedOut = bailIfBottomedOut->count > 0;

  BF_ASSERT((opts->freqTreeDepth != BF_SIZE_BAD_VALUE) ^ (opts->freqTreeOffset != BF_SIZE_BAD_VALUE));

cleanup:
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));

  return success;
}

void printNode(BfTree const *tree, BfTreeNode const *treeNode, FILE *fp) {
  (void)tree;

  BfIntervalTreeNode const *intervalTreeNode = bfTreeNodeConstToIntervalTreeNodeConst(treeNode);

  fprintf(fp, "depth = %lu", bfTreeNodeGetDepth(treeNode));
  if (bfTreeNodeIsRoot(treeNode))
    fprintf(fp, ", root");
  else
    fprintf(fp, ", index = %lu", treeNode->index);
  fprintf(fp, ": [%g, %g)", intervalTreeNode->a, intervalTreeNode->b);
  fprintf(fp, " -> %lu freqs\n", bfTreeNodeGetNumPoints(treeNode));
}

int main(int argc, char *argv[]) {
  Opts opts;
  if (!parseArgs(&opts, argc, argv))
    return EXIT_FAILURE;

  bfSeed(opts.seed);
  bfSetLogLevel(opts.logLevel);

  BF_ERROR_BEGIN() {}

  BfTrimesh *trimesh = NULL;
  BfPoints3 *verts = NULL;
  BfSize numVerts = 0;
  BfTree *rowTree = NULL;
  BfMat *L = NULL, *M = NULL;
  if (strcmp(opts.objPath, "")) {
    trimesh = bfTrimeshNewFromObjFile(opts.objPath);
    HANDLE_ERROR();

    /* Compute a finite element discretization of the Laplace-Beltrami
    * operator on `trimesh` using linear finite elements. The stiffness
    * matrix is returned in L and the mass matrix is returned in M. The
    * mass matrix isn't diagonal but this isn't too important. */
    bfTrimeshGetLboFemDiscretization(trimesh, &L, &M);
    HANDLE_ERROR();

    verts = bfTrimeshGetVerts(trimesh);

    numVerts = bfTrimeshGetNumVerts(trimesh);
    printf("using linear FEM on triangle mesh from obj file %s (%lu verts)\n", opts.objPath, numVerts);
  } else if (strcmp(opts.matPath, "")) {
    /* Load mass matrix M and stiffness matrix L from binaries */
    char rowptrPath[200];
    char colindPath[200];
    char dataPath[200];
    strcpy(rowptrPath, opts.matPath);
    strcat(rowptrPath, "/L_rowptr.bin");
    strcpy(colindPath, opts.matPath);
    strcat(colindPath, "/L_colind.bin");
    strcpy(dataPath, opts.matPath);
    strcat(dataPath, "/L_data.bin");
    L = bfMatCsrRealToMat(bfMatCsrRealNewFromBinaryFiles(rowptrPath, colindPath, dataPath));
    if (bfMatGetNumRows(L) != bfMatGetNumCols(L)) {
      printf("L must have same number of rows and columns (got %lu x %lu matrix)\n",
             bfMatGetNumRows(L), bfMatGetNumCols(L));
      exit(EXIT_FAILURE);
    }

    rowptrPath[0] = '\0';
    colindPath[0] = '\0';
    dataPath[0] = '\0';
    strcpy(rowptrPath, opts.matPath);
    strcat(rowptrPath, "/M_rowptr.bin");
    strcpy(colindPath, opts.matPath);
    strcat(colindPath, "/M_colind.bin");
    strcpy(dataPath, opts.matPath);
    strcat(dataPath, "/M_data.bin");
    M = bfMatCsrRealToMat(bfMatCsrRealNewFromBinaryFiles(rowptrPath, colindPath, dataPath));
    if (bfMatGetNumRows(M) != bfMatGetNumCols(M)) {
      printf("M must have same number of rows and columns\n");
      exit(EXIT_FAILURE);
    }

    char vertsPath[200];
    strcpy(vertsPath, opts.matPath);
    strcat(vertsPath, "/nodes.bin");
    verts = bfPoints3NewFromBinaryFile(vertsPath);

    numVerts = bfPoints3GetSize(verts);
    printf("using FEM matrices loaded from binaries in directory %s (%lu verts)\n", opts.matPath, numVerts);
  }

  if (opts.useOctree) {
    BfOctree *octree = bfOctreeNew();
    HANDLE_ERROR();

    bfOctreeInit(octree, verts, NULL, /* maxLeafSize: */ 64);
    HANDLE_ERROR();

    char const *octreeBoxesPath = "octree_boxes.txt";
    bfOctreeSaveBoxesToTextFile(octree, octreeBoxesPath);
    HANDLE_ERROR();
    printf("wrote octree cells to %s\n", octreeBoxesPath);

    rowTree = bfOctreeToTree(octree);
  } else if (opts.useFiedlerTree) {
    BfFiedlerTree *fiedlerTree = bfFiedlerTreeNewFromTrimesh(
      trimesh, /* tol: */ 1e-15, /* keepNodeTrimeshes: */ false);
    printf("built Fiedler tree\n");

    rowTree = bfFiedlerTreeToTree(fiedlerTree);
  }

  BfSize rowTreeDepth = bfTreeGetMaxDepth(rowTree);
  printf("row tree has depth %lu\n", rowTreeDepth);

  /* Find the largest eigenvalue. We need this to determine the
   * interval on which we'll build the frequency tree. */
  bfToc();
  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  HANDLE_ERROR();
  printf("computed lambda_max = %g [%.1fs]\n", lamMax, bfToc());

  /* The natural frequency of each eigenvector is the square root of
   * the associated eigenvalue. */
  BfReal freqMax = sqrt(lamMax);

  /* Set the valence of the frequency tree to match the valence of the
   * row tree: */
  BfSize k = opts.useOctree ? 8 : 2;

  /* Figure out the depth of the frequency tree: */
  BfSize freqTreeDepth = BF_SIZE_BAD_VALUE;
  if (opts.freqTreeDepth != BF_SIZE_BAD_VALUE) {
    freqTreeDepth = opts.freqTreeDepth;
  } else {
    freqTreeDepth = rowTreeDepth - opts.freqTreeOffset;
  }
  if (freqTreeDepth > rowTreeDepth) {
    printf("frequency tree depth (%lu) exceeds row tree depth\n", freqTreeDepth);
    exit(EXIT_FAILURE);
  }

  printf("building frequency tree with depth %lu (k = %lu)\n", freqTreeDepth, k);

  /* Set up the frequency tree. Note: we build the tree on the
   * frequency scale as opposed to the eigenvalue scale to preserve
   * the time-frequency product in the butterfly factorization. */
  BfIntervalTree *freqTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(freqTree, 0, freqMax, k, freqTreeDepth);
  HANDLE_ERROR();

  /* Upcast frequency tree to get the column tree */
  BfTree *colTree = bfIntervalTreeToTree(freqTree);

  FILE *fp = fopen("freqTree.txt", "w");
  bfTreeMapConst(colTree, NULL, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfTreeMapConstFunc)printNode, fp);
  fclose(fp);

  BfPoints1 *freqs = bfPoints1New();
  HANDLE_ERROR();

  bfPoints1InitEmpty(freqs, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  BfFacSpec spec = {
    .rowTree = rowTree,
    .colTree = colTree,
    .rowTreeInitDepth = opts.rowTreeOffset,
    .colTreeInitDepth = freqTreeDepth, // TODO: this is unused!
    .tol = opts.tol,
    .minNumRows = 20,
    .minNumCols = 20,
    .compareRelativeErrors = opts.compareRelativeErrors,
    .bailIfBottomedOut = opts.bailIfBottomedOut,
  };

  /* Set up the depth-first butterfly factorization streamer. We'll
   * use this below to construct the butterfly factorization
   * incrementally. */
  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, &spec);
  HANDLE_ERROR();

  if (opts.numLeafNodes != BF_SIZE_BAD_VALUE)
    printf("streaming %lu leaf nodes\n", opts.numLeafNodes);
  else if (opts.bailIfBottomedOut)
    printf("streaming until we bottom out in the space tree\n");
  else
    printf("streaming the entire eigenvector matrix (PROBABLY A BAD IDEA!)\n");

  BfReal t_total = 0;
  BfReal t_eigs = 0;

  BfFacSpan *facSpan = NULL;
  BfMat *mat = NULL;

  BfSize numStreamed = 0;
  while (!bfFacStreamerIsDone(facStreamer)) {
    BfReal t0_column = bfTime();

    BfLboFeedResult result = bfLboFeedFacStreamerNextEigenband(facStreamer, freqs, L, M);

    BfSize numBytesUncompressed = sizeof(BfReal)*numVerts*freqs->size;
    printf("- streamed %lu eigs (%1.1f%% of total) in %1.2fs\n",
           freqs->size,
           (100.0*freqs->size)/numVerts,
           result.eigenbandTime);
    printf("- uncompressed size: %.3f MB\n", numBytesUncompressed/pow(1024, 2));

    BfReal t1_column = bfTime();
    printf("- total time: %1.2fs\n", t1_column - t0_column);

    t_total += t1_column - t0_column;
    t_eigs += result.eigenbandTime;

    if (!result.success) {
      printf("* bottomed out!\n");
      break;
    }

    facSpan = bfFacStreamerGetFacSpan(facStreamer);
    mat = bfFacSpanGetMat(facSpan, BF_POLICY_VIEW);
    BfSize numBytesCompressed = bfMatNumBytes(mat);
    bfMatDelete(&mat);
    bfFacSpanDelete(&facSpan);

    printf("- compressed size:   %.3f MB\n", numBytesCompressed/pow(1024, 2));
    printf("- compression rate:  %.3f\n", (double)numBytesUncompressed/numBytesCompressed);

    HANDLE_ERROR();

    if (++numStreamed >= opts.numLeafNodes) break;
  }

  printf("finished streaming butterfly factorization\n");
  printf("- total time: %1.2fs\n", t_total);
  printf("- eigs time: %1.2fs (%0.1f%%)\n", t_eigs, 100.0*t_eigs/t_total);

  BF_ERROR_END() {}

  /** Clean up: */

  bfMatDelete(&mat);
  bfFacSpanDelete(&facSpan);
  bfFacStreamerDelete(&facStreamer);
  bfPoints1Delete(&freqs);
  bfTreeDelete(&colTree);
  bfMatDelete(&M);
  bfMatDelete(&L);
  bfTreeDelete(&rowTree);
  bfTrimeshDeinitAndDealloc(&trimesh);
}
