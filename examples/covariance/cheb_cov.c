#include <bf/assert.h>
#include <bf/cheb.h>
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

static BfReal kappa = BF_NAN;
static BfReal nu = BF_NAN;

static BfReal gamma_(BfReal lambda) {
  if (nu == 0) {
    // squared exponential spectral density function
    return exp(-kappa*lambda*lambda);
  } else {
    // Matern spectral density function, normalized so g(0) = 1
    return pow(fabs(1 + kappa*kappa*lambda), -nu/4 - 1./2);
  }
}

static BfVec *applyS(BfMat const *L, BfMat *Mfact, BfVec *w, int lumping) {
  BfVec *y    = NULL;
  BfVec *tmp1 = NULL;
  BfVec *tmp2 = NULL;
  if (lumping) {
    // Mfact is inverse diagonal sqrt C^2 = M^{-1}
    // apply C * L * C
    tmp1 = bfMatMulVec(Mfact, w);
    tmp2 = bfMatMulVec(L, tmp1);
    y    = bfMatMulVec(Mfact, tmp2);
  } else {
    // Mfact is cholesky factor C*C^T = M.
    // we want to apply C^{-1} * L * C^{-T}
    tmp1 = bfCholCsrRealSolveVec(Mfact, w); // TODO: needs to be C^{-T} not C^{-1} !
    tmp2 = bfMatMulVec(L, tmp1);
    y    = bfCholCsrRealSolveVec(Mfact, tmp2);
  }

  if (tmp1 != NULL) bfVecDelete(&tmp1);
  if (tmp2 != NULL) bfVecDelete(&tmp2);

  return y;
}

static BfVec *chebmul(BfCheb const *cheb, BfMat const *L, BfMat *Mfact, BfVec *w, int lumping) {
  BfReal const *c = cheb->c;

  BF_ASSERT(cheb->a == 0);
  BfReal lamMax = cheb->b;

  BfSize n = bfMatGetNumRows(L);
  BfVec *x = bfVecRealToVec(bfVecRealNewWithValue(n, 0));

  BfVec *y2 = bfVecCopy(w);

  bfVecDaxpy(x, c[0], y2);

  BfVec *y1 = applyS(L, Mfact, w, lumping);
  
  bfVecDscal(y1, 2/lamMax);
  bfVecDaxpy(y1, -1, w);
  bfVecDaxpy(x, c[1], y1);

  for (BfSize k = 2; k < cheb->order; ++k) {
    BfVec *y = applyS(L, Mfact, y1, lumping);
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

static BfVec *sample_z(BfCheb const *cheb, BfMat const *L, BfMat const *Mfact, int lumping) {
  BfSize n = bfMatGetNumRows(L);
  BfVec *w = bfVecRealToVec(bfVecRealNewRandn(n));
  BfVec *x = chebmul(cheb, L, Mfact, w, lumping);
  BfVec *z = NULL; 
  if (lumping) {
    z = bfMatMulVec(Mfact, x);
  } else {
    z = bfCholCsrRealSolveVec(Mfact, x); // TODO: needs to be C^{-T} not C^{-1} !
  }
  
  bfVecDelete(&w);
  bfVecDelete(&x);
  return z;
}

static BfVec *cov_matvec(BfVec *v, BfCheb const *cheb, BfMat const *L, BfMat const *Mfact, int lumping) {
  BfVec *tmp = NULL;
  if (lumping) {
    tmp = bfMatMulVec(Mfact, v);
  } else {
    tmp = bfCholCsrRealSolveVec(Mfact, v);
  }
  BfVec *tmp1 = chebmul(cheb, L, Mfact, tmp, lumping);
  bfVecDelete(&tmp);
  tmp = chebmul(cheb, L, Mfact, tmp1, lumping);
  bfVecDelete(&tmp1);
  tmp1 = NULL;
  if (lumping) {
    tmp1 = bfMatMulVec(Mfact, tmp);
  } else {
    tmp1 = bfCholCsrRealSolveVec(Mfact, tmp); // TODO: needs to be C^{-T} not C^{-1} !
  }
  bfVecDelete(&tmp);
  return tmp1;
}

static void get_factM(BfMat *M, BfMat **MfactHandle, int lumping) {
  if (lumping) {
    BfSize n = bfMatGetNumRows(M);

    BfMatCsrReal const *MCsr = bfMatConstToMatCsrRealConst(M);

    BfMatDiagReal *MLumpSqrtInv = bfMatDiagRealNew();
    bfMatDiagRealInit(MLumpSqrtInv, n, n);

    for (BfSize i = 0; i < n; ++i) {
      MLumpSqrtInv->data[i] = 0;
      for (BfSize j = MCsr->rowptr[i]; j < MCsr->rowptr[i + 1]; ++j)
        MLumpSqrtInv->data[i] += MCsr->data[j];
      MLumpSqrtInv->data[i] = 1/sqrt(MLumpSqrtInv->data[i]);
    }

    *MfactHandle = bfMatDiagRealToMat(MLumpSqrtInv);
  } else {
    BfCholCsrReal *MSqrt = bfCholCsrRealNew();
    bfCholCsrRealInit(MSqrt, M);

    *MfactHandle = MSqrt;
  }
}

int main(int argc, char const *argv[]) {
  if (argc < 5) {
    printf("usage: %s mesh.obj kappa nu num_samples [p] [lumping] \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BfSize p = argc > 5 ? strtod(argv[5], NULL) : 128;
  int lumping = argc > 6 ? atoi(argv[6]) : 0;
  kappa = atof(argv[2]);
  nu = atof(argv[3]);
  BfSize numSamples = atoi(argv[4]);

  bfSeed(0); // must seed before using PRNG

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(argv[1]);
  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  printf("triangle mesh with %lu verts\n", numVerts);

  bfToc();

  BfMat *L = NULL, *M = NULL;
  bfTrimeshGetLboFemDiscretization(trimesh, &L, &M);
  printf("set up FEM discretization [%0.1fs]\n", bfToc());

  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  printf("maximum eigenvalue: lambda = %g [%0.1fs]\n", lamMax, bfToc());

  BfCheb gammaCheb;
  bfChebInitWithDegree(&gammaCheb, p, 0, lamMax);
  bfChebInterp(&gammaCheb, gamma_, NULL);

  // /** Evaluate gamma and our Chebyshev approximation of gamma on a
  //  ** grid, and write the grid, the true gamma values, and the
  //  ** polynomial values to disk for plotting. */

  // BfSize N = 1000;
  // FILE *fp = NULL;

  char filename[50];
  // sprintf(filename, "lambda_p%i.bin", (int)p);
  // fp = fopen(filename, "w");
  // for (BfSize i = 0; i <= N; ++i) {
  //   BfReal lam = (i*lamMax)/N;
  //   fwrite(&lam, sizeof(BfReal), 1, fp);
  // }
  // fclose(fp);

  // sprintf(filename, "gamma_p%i.bin", (int)p);
  // fp = fopen(filename, "w");
  // for (BfSize i = 0; i <= N; ++i) {
  //   BfReal y = gamma_((i*lamMax)/N);
  //   fwrite(&y, sizeof(BfReal), 1, fp);
  // }
  // fclose(fp);

  // sprintf(filename, "gamma_cheb_p%i.bin", (int)p);
  // fp = fopen(filename, "w");
  // for (BfSize i = 0; i <= N; ++i) {
  //   BfReal y = bfChebEval(&gammaCheb, (i*lamMax)/N);
  //   fwrite(&y, sizeof(BfReal), 1, fp);
  // }
  // fclose(fp);

  /** Get S = sqrt(M)^-T*L*sqrt(M)^-1 */

  BfMat *Mfact = NULL;
  get_factM(M, &Mfact, lumping);

  /** Sample z once and write it out to disk for plotting. */

  BfVec *z = sample_z(&gammaCheb, L, Mfact, lumping);
  
  sprintf(filename, "z_cheb_p%lu_kappa%.1e_nu%.1e.bin", p, kappa, nu);
  bfVecSave(z, filename);
  bfVecDelete(&z);

  /** Time how long it takes to sample z numSamples times. */

  printf("computing %i samples\n", numSamples);
  bfToc();
  for (BfSize _ = 0; _ < numSamples; ++_) {
    z = sample_z(&gammaCheb, L, Mfact, lumping);
    bfVecDelete(&z);
  }
  double sampling_time = bfToc();
  printf("drew %lu samples [%0.1fs]\n", numSamples, sampling_time);

  // save factorization time and memory sizes to file
  char line[100];
  FILE *fptr;
  sprintf(filename, "performance_kappa%.1e_nu%.1e.txt", kappa, nu);
  sprintf(line, "%lu\t%.8e\n", p, sampling_time/numSamples);
  fptr = fopen(filename, "a");
  fprintf(fptr, line);
  fclose(fptr);

  /** Evaluate the covariance function with respect to a fixed point
   ** on the mesh. */

  BfVec *e = bfVecRealToVec(bfVecRealNewStdBasis(numVerts, 0));
  BfVec *c = cov_matvec(e, &gammaCheb, L, Mfact, lumping);
  sprintf(filename, "c_cheb_p%lu_kappa%.1e_nu%.1e.bin", p, kappa, nu);
  bfVecSave(c, filename);

  /* Compute and store covariance matrix vector products */

  bfSeed(0);
  BfSize s = numSamples;
  BfMatDenseReal *randvecs = bfMatDenseRealNewZeros(numVerts, s);
  BfMatDenseReal *matvecs = bfMatDenseRealNewZeros(numVerts, s);

  printf("computing %lu matvecs with covariance\n", s);
  for (BfSize j = 0; j < s; ++j) {
      BfVec *x = bfVecRealToVec(bfVecRealNewRandn(numVerts));

      // Set column of inputs:
      bfMatDenseRealSetCol(randvecs, j, x);

      // apply covariance matrix
      BfVec *tmp1 = cov_matvec(x, &gammaCheb, L, Mfact, lumping);

      // Set column of results:
      bfMatDenseRealSetCol(matvecs, j, tmp1);

      // Clean up:
      bfVecDelete(&x);
      bfVecDelete(&tmp1);
  }

  sprintf(filename, "randvecs_cheb_p%i_kappa%.1e_nu%.1e.bin", p, kappa, nu);
  bfMatDenseRealSave(randvecs, filename);

  sprintf(filename, "matvecs_cheb_p%i_kappa%.1e_nu%.1e.bin", p, kappa, nu);
  bfMatDenseRealSave(matvecs, filename);

  printf("cleaning up\n");

  /** Clean up: */

  bfTrimeshDeinitAndDealloc(&trimesh);
  bfMatDelete(&L);
  bfMatDelete(&M);
  bfChebDeinit(&gammaCheb);
  bfVecDelete(&c);
  bfMatDelete(&Mfact);
}
