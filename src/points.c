#include <bf/points.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <bf/bbox.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_real.h>
#include <bf/vec_real.h>

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q) {
  return hypot(q[0] - p[0], q[1] - p[1]);
}

BfReal bfPoint3Dist(BfPoint3 const p, BfPoint3 const q) {
  BfVector3 u = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

void bfPoint3Sub(BfPoint3 const u, BfPoint3 const v, BfVector3 w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
}

void bfPoint3GetPointOnRay(BfPoint3 const r0, BfVector3 const dr, BfReal t, BfPoint3 rt) {
  rt[0] = r0[0] + t*dr[0];
  rt[1] = r0[1] + t*dr[1];
  rt[2] = r0[2] + t*dr[2];
}

void bfPoint3Copy(BfPoint3 x, BfPoint3 const y) {
  x[0] = y[0];
  x[1] = y[1];
  x[2] = y[2];
}

BfPoints2 const *bfPoints2ConstViewFromMat(BfMat const *mat) {
  return bfPoints2ConstViewFromMatDenseReal(bfMatConstToMatDenseRealConst(mat));
}

BfPoints2 const *bfPoints2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal) {
  BEGIN_ERROR_HANDLING();

  BfPoints2 *points = NULL;

  BfMat const *mat = bfMatDenseRealConstToMatConst(matDenseReal);
  BfMatDense const *matDense = bfMatDenseRealConstToMatDenseConst(matDenseReal);

  if (bfMatDenseRealGetNumCols(mat) != 2)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfMatDenseGetColStride(matDense) != 1)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  points = malloc(sizeof(BfPoints2));
  points->data = (BfPoint2 *)matDenseReal->data;
  points->size = bfMatDenseRealGetNumRows(mat);

  END_ERROR_HANDLING() {
    free(points);
    points = NULL;
  }

  return points;
}

/** Implementation: Points1 */

BfPoints1 *bfPoints1New() {
  BEGIN_ERROR_HANDLING();

  BfPoints1 *points = malloc(sizeof(BfPoints1));
  if (points == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  points->data = NULL;
  points->size = BF_SIZE_BAD_VALUE;
  points->capacity = BF_SIZE_BAD_VALUE;

  END_ERROR_HANDLING() {
    points = NULL;
  }

  return points;
}

void bfPoints1InitEmpty(BfPoints1 *points, BfSize capacity) {
  BEGIN_ERROR_HANDLING();

  points->data = malloc(capacity*sizeof(BfPoint1));
  if (points->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  points->size = 0;
  points->capacity = capacity;
  points->isView = false;

  END_ERROR_HANDLING() {
    bfPoints1Deinit(points);
  }
}

void bfPoints1InitViewFromVecReal(BfPoints1 *points, BfVecReal const *vecReal) {
  BfVec const *vec = bfVecRealConstToVecConst(vecReal);

  points->data = vecReal->data;
  points->size = vec->size;
  points->capacity = BF_SIZE_BAD_VALUE;
  points->isView = true;
}

void bfPoints1Deinit(BfPoints1 *points) {
  if (!points->isView)
    free(points->data);
  points->data = NULL;
  points->size = BF_SIZE_BAD_VALUE;
  points->capacity = BF_SIZE_BAD_VALUE;
}

bool bfPoints1IsSorted(BfPoints1 const *points) {
  for (BfSize i = 1; i < points->size; ++i)
    if (points->data[i - 1] > points->data[i])
      return false;
  return true;
}

void bfPoints1InsertPointsSorted(BfPoints1 *points, BfPoints1 const *newPoints, BfSize **posPtr) {
  BEGIN_ERROR_HANDLING();

  if (!bfPoints1IsSorted(points))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!bfPoints1IsSorted(newPoints))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (posPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*posPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize *pos = malloc(newPoints->size*sizeof(BfSize));
  if (pos == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize newSize = points->size + newPoints->size;

  BfReal *newData = malloc(newSize*sizeof(BfReal));
  if (newData == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize i = 0, j = 0, k = 0;
  while ((i < points->size || j < newPoints->size) && k < newSize) {
    if (i < points->size && j < newPoints->size) {
      if (points->data[i] < newPoints->data[j]) {
        newData[k++] = points->data[i++];
      } else {
        pos[j] = k;
        newData[k++] = newPoints->data[j++];
      }
    } else if (i < points->size) {
      newData[k++] = points->data[i++];
    } else if (j < newPoints->size) {
      pos[j] = k;
      newData[k++] = newPoints->data[j++];
    } else {
      assert(false);
    }
  }
  assert(i == points->size && j == newPoints->size && k == newSize);

  if (!points->isView)
    free(points->data);
  points->data = newData;
  points->size = newSize;
  points->capacity = newSize;

  END_ERROR_HANDLING() {
    free(pos);
    free(newData);
    *posPtr = NULL;
  }

  *posPtr = pos;
}

/** Implementation: Points2 */

BfPoints2 bfGetUninitializedPoints2() {
  return (BfPoints2) {.data = NULL, .size = 0};
}

void bfInitEmptyPoints2(BfPoints2 *points, BfSize numPoints) {
  if (numPoints == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  points->size = numPoints;

  points->data = malloc(numPoints*sizeof(BfPoint2));
  if (points->data == NULL)
    bfSetError(BF_ERROR_MEMORY_ERROR);
}

void bfReadPoints2FromFile(char const *path, BfPoints2 *points) {
  BEGIN_ERROR_HANDLING();

  /* open the file for reading */
  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the size of the file */
  fseek(fp, 0, SEEK_END);
  BfSize size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  /* make sure the binary file is the right size */
  if (size % sizeof(BfPoint2) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the number of points */
  points->size = size/sizeof(BfPoint2);

  /* allocate space for the points */
  points->data = malloc(size);
  if (points->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* read them in */
  fread(points->data, sizeof(BfPoint2), points->size, fp);
  // TODO: error-handling

  END_ERROR_HANDLING() {
    free(points->data);
  }

  fclose(fp);
}

void bfFreePoints2(BfPoints2 *points) {
  free(points->data);
}

bool bfPoints2Initialized(BfPoints2 const *points) {
  return points->data != NULL && points->size != 0;
}

BfBbox2 bfPoints2GetBoundingBox(BfPoints2 const *points) {
  BfBbox2 bbox = bfGetEmptyBbox2();

  BfSize numPoints = points->size;

  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal const *point = points->data[i];

    BfReal x = point[0], y = point[1];

    bbox.min[0] = fmin(bbox.min[0], x);
    bbox.max[0] = fmax(bbox.max[0], x);
    bbox.min[1] = fmin(bbox.min[1], y);
    bbox.max[1] = fmax(bbox.max[1], y);
  }

  return bbox;
}

void bfGetPointsByIndex(BfPoints2 const *points,
                        BfSize numInds, BfSize const *inds,
                        BfPoints2 *indexedPoints)
{
  BEGIN_ERROR_HANDLING();

  bfInitEmptyPoints2(indexedPoints, numInds);
  HANDLE_ERROR();

  BfPoint2 const *point = (BfPoint2 const *)points->data;
  BfPoint2 *indexedPoint = indexedPoints->data;

  for (BfSize i = 0, j; i < numInds; ++i) {
    j = inds[i];
    indexedPoint[i][0] = point[j][0];
    indexedPoint[i][1] = point[j][1];
  }

  END_ERROR_HANDLING() {
    bfFreePoints2(indexedPoints);
  }
}

void bfPrintPoints2(BfPoints2 const *points) {
  for (BfSize i = 0; i < points->size; ++i)
    printf("points[%lu] = (%g, %g)\n",i,points->data[i][0],points->data[i][1]);
}

void bfSavePoints2(BfPoints2 const *points, char const *path) {
  BEGIN_ERROR_HANDLING();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(points->data, points->size, sizeof(BfPoint2), fp);
  // TODO: error-handling

  END_ERROR_HANDLING() {}

  fclose(fp);
}

/* Compute all pairwise distances between the point sets `X` and
 * `Y`. If `|X| == m` and `|Y| == n`, then the result is a length
 * `m*n` array of `BfReal`s where the `m*i + j` element is the
 * distance between `X[i]` and `Y[j]`.
 *
 * It is the responsibility of the caller to free the returned
 * array. */
BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y) {
  BEGIN_ERROR_HANDLING();

  BfSize m = X->size, n = Y->size;

  BfReal *r = malloc(m*n*sizeof(BfReal));
  if (r == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    BfReal const *x = &X->data[i][0];
    for (BfSize j = 0; j < n; ++j) {
      BfReal const *y = &Y->data[j][0];
      r[k++] = hypot(y[0] - x[0], y[1] - x[1]);
    }
  }

  END_ERROR_HANDLING() {}

  return r;
}

/** Implementation: Points3 */

void bfPoints3InitEmpty(BfPoints3 *points, BfSize numPoints) {
  if (numPoints == 0) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  points->size = numPoints;

  points->data = malloc(numPoints*sizeof(BfPoint3));
  if (points->data == NULL)
    bfSetError(BF_ERROR_MEMORY_ERROR);
}

void bfPoints3InitFromBinaryFile(BfPoints3 *points, char const *path) {
  BEGIN_ERROR_HANDLING();

  /* open the file for reading */
  FILE *fp = fopen(path, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the size of the file */
  fseek(fp, 0, SEEK_END);
  BfSize size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  /* make sure the binary file is the right size */
  if (size % sizeof(BfPoint3) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* get the number of points */
  points->size = size/sizeof(BfPoint3);

  /* allocate space for the points */
  points->data = malloc(size);
  if (points->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* read them in */
  fread(points->data, sizeof(BfPoint3), points->size, fp);
  // TODO: error-handling

  END_ERROR_HANDLING() {
    free(points->data);
  }

  fclose(fp);
}

void bfPoints3Deinit(BfPoints3 *points) {
  (void)points;
  assert(false);
}

BfBoundingBox3 bfPoints3GetBoundingBox(BfPoints3 const *points) {
  BfBoundingBox3 boundingBox;
  bfBoundingBox3InitEmpty(&boundingBox);

  BfSize numPoints = points->size;

  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal const *point = points->data[i];
    for (BfSize j = 0; j < 3; ++j) {
      boundingBox.min[j] = fmin(boundingBox.min[j], point[j]);
      boundingBox.max[j] = fmax(boundingBox.max[j], point[j]);
    }
  }

  return boundingBox;
}

void bfPoints3GetByIndex(BfPoints3 const *points, BfSize numInds, BfSize const *inds, BfPoints3 *indexedPoints) {
  BEGIN_ERROR_HANDLING();

  bfPoints3InitEmpty(indexedPoints, numInds);
  HANDLE_ERROR();

  for (BfSize i = 0; i < numInds; ++i)
    bfPoint3Copy(indexedPoints->data[i], points->data[inds[i]]);

  END_ERROR_HANDLING() {
    bfPoints3Deinit(indexedPoints);
  }
}
