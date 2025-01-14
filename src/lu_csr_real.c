#include <bf/lu_csr_real.h>

#include <limits.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_csr_real.h>
#include <bf/mem.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#if defined BF_LINUX
#  include <suitesparse/umfpack.h>
#elif defined BF_DARWIN
#  include <umfpack.h>
#else
#  error "Unknown platform"
#endif

/** Interface: Lu */

static BfLuVtable LU_VTABLE = {
  .Solve = (__typeof__(&bfLuSolve))bfLuCsrRealSolve,
  .SolveLower = (__typeof__(&bfLuSolveLower))bfLuCsrRealSolveLower,
  .SolveUpper = (__typeof__(&bfLuSolveUpper))bfLuCsrRealSolveUpper,
  .Scale = (__typeof__(&bfLuScale))bfLuCsrRealScale,
  .SolveVec = (__typeof__(&bfLuSolveVec))bfLuCsrRealSolveVec,
  .SolveLowerVec = (__typeof__(&bfLuSolveLowerVec))bfLuCsrRealSolveLowerVec,
  .SolveUpperVec = (__typeof__(&bfLuSolveUpperVec))bfLuCsrRealSolveUpperVec,
  .ScaleVec = (__typeof__(&bfLuScaleVec))bfLuCsrRealScaleVec,
};

BfMat *bfLuCsrRealSolve(BfLuCsrReal const *luCsrReal, BfMat const *B) {
  (void)luCsrReal;
  (void)B;
  BF_DIE();
}

BfMat *bfLuCsrRealSolveLower(BfLuCsrReal const *luCsrReal, BfMat const *B, bool permute) {
  (void)luCsrReal;
  (void)B;
  (void)permute;
  BF_DIE();
}

BfMat *bfLuCsrRealSolveUpper(BfLuCsrReal const *luCsrReal, BfMat const *B, bool permute) {
  (void)luCsrReal;
  (void)B;
  (void)permute;
  BF_DIE();
}

BfMat *bfLuCsrRealScale(BfLuCsrReal const *luCsrReal, BfMat const *B) {
  (void)luCsrReal;
  (void)B;
  BF_DIE();
}

static BfVec *solve_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b) {
  BF_ERROR_BEGIN();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    UMFPACK_At,
    luCsrReal->rowptr,
    luCsrReal->colind,
    luCsrReal->data,
    x->data,
    b->data,
    luCsrReal->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ERROR_END() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealSolveVec(BfLuCsrReal const *luCsrReal, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solve_vecReal(luCsrReal, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *solveLowerVec_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b, bool permute) {
  BF_ERROR_BEGIN();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    permute ? UMFPACK_Q_Ut : UMFPACK_Ut,
    luCsrReal->rowptr,
    luCsrReal->colind,
    luCsrReal->data,
    x->data,
    b->data,
    luCsrReal->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ERROR_END() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealSolveLowerVec(BfLuCsrReal const *luCsrReal, BfVec const *b, bool permute) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solveLowerVec_vecReal(luCsrReal, bfVecConstToVecRealConst(b), permute);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *solveUpperVec_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b, bool permute) {
  BF_ERROR_BEGIN();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_solve(
    permute ? UMFPACK_Lt_P : UMFPACK_Lt,
    luCsrReal->rowptr,
    luCsrReal->colind,
    luCsrReal->data,
    x->data,
    b->data,
    luCsrReal->numeric,
    NULL,
    NULL);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ERROR_END() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealSolveUpperVec(BfLuCsrReal const *luCsrReal, BfVec const *b, bool permute) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return solveUpperVec_vecReal(luCsrReal, bfVecConstToVecRealConst(b), permute);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfVec *scale_vecReal(BfLuCsrReal const *luCsrReal, BfVecReal const *b) {
  BF_ERROR_BEGIN();

  BfVecReal *x = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(x, b->super.size);
  HANDLE_ERROR();

  int status = umfpack_di_scale(x->data, b->data, luCsrReal->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ERROR_END() {
    bfVecRealDeinitAndDealloc(&x);
  }

  return bfVecRealToVec(x);
}

BfVec *bfLuCsrRealScaleVec(BfLuCsrReal const *luCsrReal, BfVec const *b) {
  switch (bfVecGetType(b)) {
  case BF_TYPE_VEC_REAL:
    return scale_vecReal(luCsrReal, bfVecConstToVecRealConst(b));
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

/** Upcasting: LuCsrReal -> Lu */

/** Downcasting: Lu -> LuCsrReal */

/** Implementation: LuCsrReal */

BfLuCsrReal *bfLuCsrRealNew() {
  BF_ERROR_BEGIN();

  BfLuCsrReal *luCsrReal = bfMemAlloc(1, sizeof(BfLuCsrReal));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return luCsrReal;
}

static void checkUmfpackStatus(int status) {
  switch (status) {
  case UMFPACK_OK:
    break;
  case UMFPACK_ERROR_n_nonpositive:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  case UMFPACK_ERROR_invalid_matrix:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  case UMFPACK_ERROR_out_of_memory:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  case UMFPACK_ERROR_argument_missing:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  case UMFPACK_ERROR_internal_error:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  default:
    BF_DIE();
  }
}

void bfLuCsrRealInit(BfLuCsrReal *luCsrReal, BfMat const *mat) {
  BF_ERROR_BEGIN();

  bfLuInit(&luCsrReal->super, &LU_VTABLE);

  BfMatCsrReal const *matCsrReal = bfMatConstToMatCsrRealConst(mat);
  HANDLE_ERROR();

  BfSize n = bfMatGetNumRows(mat);
  if (n != bfMatGetNumCols(mat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Can't exceed range of int since UMFPACK uses ints for indices. */
  if (n > INT_MAX)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Each range of column indices must be sorted in ascending order: */
  for (BfSize i = 0; i < n; ++i) {
    BfSize k0 = matCsrReal->rowptr[i];
    BfSize k1 = matCsrReal->rowptr[i + 1];
    if (!bfSizeIsStrictlyIncreasing(k1 - k0, &matCsrReal->colind[k0]))
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  BfSize nnz = matCsrReal->rowptr[n];

  luCsrReal->rowptr = bfMemAlloc(n + 1, sizeof(int));
  if (luCsrReal->rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i <= n; ++i)
    luCsrReal->rowptr[i] = matCsrReal->rowptr[i];

  luCsrReal->colind = bfMemAlloc(nnz, sizeof(int));
  if (luCsrReal->colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < nnz; ++i)
    luCsrReal->colind[i] = matCsrReal->colind[i];

  luCsrReal->data = bfMemAlloc(nnz, sizeof(BfReal));
  if (luCsrReal->data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfMemCopy(matCsrReal->data, nnz, sizeof(BfReal), luCsrReal->data);

  /* Use UMFPACK to compute the LU decomposition of `mat`. Since
   * UMFPACK operates on CSC sparse matrices, and since `mat` is a CSR
   * sparse matrix, we actually compute the LU decomposition A^T =
   * U^T*L*T. */

  int status;

  status = umfpack_di_symbolic(
    n, n, luCsrReal->rowptr, luCsrReal->colind, luCsrReal->data,
    &luCsrReal->symbolic, NULL, NULL);
  checkUmfpackStatus(status);

  status = umfpack_di_numeric(
    luCsrReal->rowptr, luCsrReal->colind, luCsrReal->data,
    luCsrReal->symbolic, &luCsrReal->numeric, NULL, NULL);
  checkUmfpackStatus(status);

  BF_ERROR_END() {
    bfLuCsrRealDeinit(luCsrReal);
  }
}

void bfLuCsrRealDeinit(BfLuCsrReal *luCsrReal) {
  bfMemFree(luCsrReal->rowptr);
  bfMemFree(luCsrReal->colind);
  bfMemFree(luCsrReal->data);

  luCsrReal->rowptr = NULL;
  luCsrReal->colind = NULL;
  luCsrReal->data = NULL;

  umfpack_di_free_numeric(&luCsrReal->numeric);
  luCsrReal->numeric = NULL;

  umfpack_di_free_symbolic(&luCsrReal->symbolic);
  luCsrReal->symbolic = NULL;
}

void bfLuCsrRealDealloc(BfLuCsrReal **luCsrReal) {
  bfMemFree(*luCsrReal);
  *luCsrReal = NULL;
}

void bfLuCsrRealDeinitAndDealloc(BfLuCsrReal **luCsrReal) {
  bfLuCsrRealDeinit(*luCsrReal);
  bfLuCsrRealDealloc(luCsrReal);
}
void bfLuCsrRealDump(BfLuCsrReal const *luCsrReal) {
  BF_ERROR_BEGIN();

  int *Lp = NULL;
  int *Lj = NULL;
  double *Lx = NULL;
  int *Up = NULL;
  int *Ui = NULL;
  double *Ux = NULL;
  int *P = NULL;
  int *Q = NULL;
  double *Dx = NULL;
  double *Rs = NULL;

  int status = 0;

  int lnz, unz, n_row, n_col, nz_udiag;
  status = umfpack_di_get_lunz(&lnz, &unz, &n_row, &n_col, &nz_udiag, luCsrReal->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BF_ASSERT(n_row == n_col);
  int n = n_row;

  Lp = bfMemAlloc(n + 1, sizeof(int));
  if (Lp == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Lj = bfMemAlloc(lnz, sizeof(int));
  if (Lj == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Lx = bfMemAlloc(lnz, sizeof(double));
  if (Lx == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Up = bfMemAlloc(n + 1, sizeof(int));
  if (Up == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Ui = bfMemAlloc(unz, sizeof(int));
  if (Ui == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Ux = bfMemAlloc(unz, sizeof(double));
  if (Ux == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  P = bfMemAlloc(n, sizeof(int));
  if (P == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Q = bfMemAlloc(n, sizeof(int));
  if (Q == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Dx = bfMemAlloc(n, sizeof(double));
  if (Dx == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  Rs = bfMemAlloc(n, sizeof(double));
  if (Rs == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  int do_recip;
  status = umfpack_di_get_numeric(
    Lp, Lj, Lx, Up, Ui, Ux, P, Q, Dx, &do_recip, Rs, luCsrReal->numeric);
  if (status)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  FILE *fp = NULL;

  fp = fopen("Lp.bin", "w");
  fwrite(Lp, sizeof(int), n + 1, fp);
  fclose(fp);

  fp = fopen("Lj.bin", "w");
  fwrite(Lj, sizeof(int), lnz, fp);
  fclose(fp);

  fp = fopen("Lx.bin", "w");
  fwrite(Lx, sizeof(double), lnz, fp);
  fclose(fp);

  fp = fopen("Up.bin", "w");
  fwrite(Up, sizeof(int), n + 1, fp);
  fclose(fp);

  fp = fopen("Ui.bin", "w");
  fwrite(Ui, sizeof(int), unz, fp);
  fclose(fp);

  fp = fopen("Ux.bin", "w");
  fwrite(Ux, sizeof(double), unz, fp);
  fclose(fp);

  fp = fopen("P.bin", "w");
  fwrite(P, sizeof(int), n, fp);
  fclose(fp);

  fp = fopen("Q.bin", "w");
  fwrite(Q, sizeof(int), n, fp);
  fclose(fp);

  fp = fopen("Dx.bin", "w");
  fwrite(Dx, sizeof(double), n, fp);
  fclose(fp);

  fp = fopen("Rs.bin", "w");
  fwrite(Rs, sizeof(double), n, fp);
  fclose(fp);

  BF_ERROR_END() {}

  bfMemFree(Lp);
  bfMemFree(Lj);
  bfMemFree(Lx);
  bfMemFree(Up);
  bfMemFree(Ui);
  bfMemFree(Ux);
  bfMemFree(P);
  bfMemFree(Q);
  bfMemFree(Dx);
  bfMemFree(Rs);
}
