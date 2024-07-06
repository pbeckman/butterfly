#include <bf/assert.h>
#include <bf/cheb.h>
#include <bf/chol_csr_real.h>
#include <bf/const.h>
#include <bf/linalg.h>
#include <bf/mat_csr_real.h>
#include <bf/mat_product.h>
#include <bf/mat_dense_real.h>
#include <bf/rand.h>
#include <bf/trimesh.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct {
  BfSize p;
  bool lumping;
  BfReal kappa;
  BfReal nu;
  BfSize numSamples;
  BfTrimesh *trimesh;
  BfSize n;
  BfMat *L;
  BfMat *M;
  BfReal lamMax;
  BfCheb gammaCheb;
  union {
    BfMat *MFactLumped;
    BfChol *MFactCholesky;
  };
} context;

static BfReal gamma_(BfReal lambda) {
  if (context.nu == 0) {
    // squared exponential spectral density function
    return exp(-context.kappa*lambda*lambda);
  } else {
    // Matern spectral density function, normalized so g(0) = 1
    return pow(fabs(1 + pow(context.kappa, 2)*lambda), -context.nu/4 - 1./2);
  }
}

static BfVec *applyS(BfVec *w) {
  BfVec *y    = NULL;
  BfVec *tmp1 = NULL;
  BfVec *tmp2 = NULL;
  if (context.lumping) {
    // MFact is inverse diagonal sqrt C^2 = M^{-1}
    // apply C * L * C
    tmp1 = bfMatMulVec(context.MFactLumped, w);
    tmp2 = bfMatMulVec(context.L, tmp1);
    y    = bfMatMulVec(context.MFactLumped, tmp2);
  } else {
    // MFact is cholesky factor C*C^T = M.
    // we want to apply C^{-1} * L * C^{-T}
    tmp1 = bfCholFacSolveVec(context.MFactCholesky, w, /* transposed: */ true);
    tmp2 = bfMatMulVec(context.L, tmp1);
    y    = bfCholFacSolveVec(context.MFactCholesky, tmp2, /* transposed: */ false);
  }

  if (tmp1 != NULL) bfVecDelete(&tmp1);
  if (tmp2 != NULL) bfVecDelete(&tmp2);

  return y;
}

static BfVec *chebmul(BfVec *w) {
  BfReal const *c = context.gammaCheb.c;

  BF_ASSERT(context.gammaCheb.a == 0);
  BfReal lamMax = context.gammaCheb.b;

  BfVec *x = bfVecRealToVec(bfVecRealNewWithValue(context.n, 0));

  BfVec *y2 = bfVecCopy(w);

  bfVecDaxpy(x, c[0], y2);

  BfVec *y1 = applyS(w);

  bfVecDscal(y1, 2/lamMax);
  bfVecDaxpy(y1, -1, w);
  bfVecDaxpy(x, c[1], y1);

  for (BfSize k = 2; k < context.gammaCheb.order; ++k) {
    BfVec *y = applyS(y1);
    bfVecDscal(y, 4/lamMax);
    bfVecDaxpy(y, -2, y1);
    bfVecDaxpy(y, -1, y2);

    bfVecDaxpy(x, c[k], y);

    bfVecDelete(&y2);
    y2 = y1;
    y1 = y;
  }

  if (y1 != NULL) bfVecDelete(&y1);
  if (y2 != NULL) bfVecDelete(&y2);

  return x;
}

static BfVec *sample_z(void) {
  BfVec *w = bfVecRealToVec(bfVecRealNewRandn(context.n));
  BfVec *x = chebmul(w);
  BfVec *z = NULL;
  if (context.lumping) {
    z = bfMatMulVec(context.MFactLumped, x);
  } else {
    z = bfCholFacSolveVec(context.MFactCholesky, x, /* transposed: */ true);
  }

  bfVecDelete(&w);
  bfVecDelete(&x);
  return z;
}

static BfVec *cov_matvec(BfVec *v) {
  BfVec *tmp = NULL;
  if (context.lumping) {
    tmp = bfMatMulVec(context.MFactLumped, v);
  } else {
    tmp = bfCholFacSolveVec(context.MFactCholesky, v, /* transposed: */ false);
  }
  BfVec *tmp1 = chebmul(tmp);
  bfVecDelete(&tmp);
  tmp = chebmul(tmp1);
  bfVecDelete(&tmp1);
  tmp1 = NULL;
  if (context.lumping) {
    tmp1 = bfMatMulVec(context.MFactLumped, tmp);
  } else {
    BF_DIE();
    tmp1 = bfCholFacSolveVec(context.MFactCholesky, tmp, /* transposed: */ true);
  }
  bfVecDelete(&tmp);
  return tmp1;
}

int main(int argc, char const *argv[]) {
  bfSeed(0); // must seed before using PRNG

  if (argc < 5) {
    printf("usage: %s mesh.obj kappa nu num_samples [p] [lumping] \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  context.p = argc > 5 ? strtod(argv[5], NULL) : 128;
  context.lumping = argc > 6 ? atoi(argv[6]) : 0;
  context.kappa = atof(argv[2]);
  context.nu = atof(argv[3]);
  context.numSamples = atoi(argv[4]);
  context.trimesh = bfTrimeshNewFromObjFile(argv[1]);
  context.n = bfTrimeshGetNumVerts(context.trimesh);

  printf("triangle mesh with %lu verts\n", context.n);

  bfToc();

  bfTrimeshGetLboFemDiscretization(context.trimesh, &context.L, &context.M);
  printf("set up FEM discretization [%0.1fs]\n", bfToc());

  context.lamMax = bfGetMaxEigenvalue(context.L, context.M);
  printf("maximum eigenvalue: lambda = %g [%0.1fs]\n", context.lamMax, bfToc());

  bfChebInitWithDegree(&context.gammaCheb, context.p, 0, context.lamMax);
  bfChebInterp(&context.gammaCheb, gamma_, NULL);

  /** Initialize MFact, taking into account whether we're doing mass lumping: */

  if (context.lumping) {
    BfMatCsrReal const *MCsr = bfMatConstToMatCsrRealConst(context.M);

    BfMatDiagReal *MLumpSqrtInv = bfMatDiagRealNew();
    bfMatDiagRealInit(MLumpSqrtInv, context.n, context.n);

    for (BfSize i = 0; i < context.n; ++i) {
      MLumpSqrtInv->data[i] = 0;
      for (BfSize j = MCsr->rowptr[i]; j < MCsr->rowptr[i + 1]; ++j)
        MLumpSqrtInv->data[i] += MCsr->data[j];
      MLumpSqrtInv->data[i] = 1/sqrt(MLumpSqrtInv->data[i]);
    }

    context.MFactLumped = bfMatDiagRealToMat(MLumpSqrtInv);
  } else {
    BfCholCsrReal *MSqrt = bfCholCsrRealNew();
    bfCholCsrRealInit(MSqrt, context.M);
    context.MFactCholesky = bfCholCsrRealToChol(MSqrt);
  }

  /** Evaluate gamma and our Chebyshev approximation of gamma on a
   ** grid, and write the grid, the true gamma values, and the
   ** polynomial values to disk for plotting. */

  BfSize N = 1000;
  FILE *fp = NULL;

  char filename[50];
  sprintf(filename, "lambda_p%lu.bin", context.p);

  fp = fopen(filename, "wb");
  for (BfSize i = 0; i <= N; ++i) {
    BfReal lam = (i*context.lamMax)/N;
    fwrite(&lam, sizeof(BfReal), 1, fp);
  }
  fclose(fp);

  sprintf(filename, "gamma_p%lu.bin", context.p);
  fp = fopen(filename, "wb");
  for (BfSize i = 0; i <= N; ++i) {
    BfReal y = gamma_((i*context.lamMax)/N);
    fwrite(&y, sizeof(BfReal), 1, fp);
  }
  fclose(fp);

  sprintf(filename, "gamma_cheb_p%lu.bin", context.p);
  fp = fopen(filename, "wb");
  for (BfSize i = 0; i <= N; ++i) {
    BfReal y = bfChebEval(&context.gammaCheb, (i*context.lamMax)/N);
    fwrite(&y, sizeof(BfReal), 1, fp);
  }
  fclose(fp);

  /** Sample z once and write it out to disk for plotting. */

  BfVec *z = sample_z();

  sprintf(filename, "z_cheb_p%lu_kappa%.1e_nu%.1e.bin", context.p, context.kappa, context.nu);
  bfVecSave(z, filename);
  bfVecDelete(&z);

  /** Time how long it takes to sample z numSamples times. */

  printf("computing %lu samples\n", context.numSamples);
  bfToc();
  for (BfSize _ = 0; _ < context.numSamples; ++_) {
    z = sample_z();
    bfVecDelete(&z);
  }
  double sampling_time = bfToc();
  printf("drew %lu samples [%0.1fs]\n", context.numSamples, sampling_time);

  // save factorization time and memory sizes to file
  char line[100];
  FILE *fptr;
  sprintf(filename, "performance_kappa%.1e_nu%.1e.txt", context.kappa, context.nu);
  sprintf(line, "%lu\t%.8e\n", context.p, sampling_time/context.numSamples);
  fptr = fopen(filename, "a");
  fprintf(fptr, line);
  fclose(fptr);

  /** Evaluate the covariance function with respect to a fixed point
   ** on the mesh. */

  BfVec *e = bfVecRealToVec(bfVecRealNewStdBasis(context.n, 0));
  BfVec *c = cov_matvec(e);
  sprintf(filename, "c_cheb_p%lu_kappa%.1e_nu%.1e.bin", context.p, context.kappa, context.nu);
  bfVecSave(c, filename);

  /* Compute and store covariance matrix vector products */


  BfMatDenseReal *randvecs = bfMatDenseRealNewZeros(context.n, context.numSamples);
  BfMatDenseReal *matvecs = bfMatDenseRealNewZeros(context.n, context.numSamples);

  printf("computing %lu matvecs with covariance\n", context.numSamples);
  for (BfSize j = 0; j < context.numSamples; ++j) {
      BfVec *x = bfVecRealToVec(bfVecRealNewRandn(context.n));

      // Set column of inputs:
      bfMatDenseRealSetCol(randvecs, j, x);

      // apply covariance matrix
      BfVec *tmp1 = cov_matvec(x);

      // Set column of results:
      bfMatDenseRealSetCol(matvecs, j, tmp1);

      // Clean up:
      bfVecDelete(&x);
      bfVecDelete(&tmp1);
  }

  sprintf(filename, "randvecs_cheb_p%lu_kappa%.1e_nu%.1e.bin", context.p, context.kappa, context.nu);
  bfMatDenseRealSave(randvecs, filename);

  sprintf(filename, "matvecs_cheb_p%lu_kappa%.1e_nu%.1e.bin", context.p, context.kappa, context.nu);
  bfMatDenseRealSave(matvecs, filename);

  printf("cleaning up\n");

  /** Clean up: */

  bfTrimeshDeinitAndDealloc(&context.trimesh);
  bfMatDelete(&context.L);
  bfMatDelete(&context.M);
  bfChebDeinit(&context.gammaCheb);
  if (context.lumping)
    bfMatDelete(&context.MFactLumped);
  else
    bfCholDelete(&context.MFactCholesky);
}
