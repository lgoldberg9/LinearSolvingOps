/*
  Last modified by Logan Goldberg on 3/4/2018
 */
#include "MATRIXOPS.c"

/* C module for basic vector and array operations. */
void swapInt(uint32_t* i1, uint32_t* i2) {
  int t = *i1;
  *i1 = *i2;
  *i2 = t;
}

void swapDouble(double* d1, double* d2) {
  double t = *d1;
  *d1 = *d2;
  *d2 = t;
}

enum MatErrType CopyMat(Matrix MatB, Matrix MatA,
			uint32_t numRows, uint32_t numCols) {
  
  if (!(checkRowCol(MatA, numRows, numCols) &&
	checkRowCol(MatB, numRows, numCols))) {
    return matErr_Size;
  }

  for (uint32_t row = 0; row < numRows; row++) {
    for (uint32_t col = 0; col < numCols; col++) {
      MAT(MatB, row, col) = MAT(MatA, row, col);
    }
  }
  return matErr_None;
}

enum MatErrType AddMat(Matrix MatC, Matrix MatA, Matrix MatB,
		       uint32_t numRows, uint32_t numCols) {
  
  if (!(checkRowCol(MatA, numRows, numCols) &&
	checkRowCol(MatB, numRows, numCols) &&
	checkRowCol(MatC, numRows, numCols))) {
    return matErr_Size;
  }

  
  for (uint32_t row = 0; row < numRows; row++) {
    for (uint32_t col = 0; col < numCols; col++) {
      MAT(MatC, row, col) = MAT(MatA, row, col) + MAT(MatB, row, col);
    }
  }
  return matErr_None;
}

enum MatErrType SubMat(Matrix MatC, Matrix MatA, Matrix MatB,
		       uint32_t numRows, uint32_t numCols) {
  
  if (!(checkRowCol(MatA, numRows, numCols) &&
	checkRowCol(MatB, numRows, numCols) &&
	checkRowCol(MatC, numRows, numCols))) {
    return matErr_Size;
  }

  for (uint32_t row = 0; row < numRows; row++) {
    for (uint32_t col = 0; col < numCols; col++) {
      MAT(MatC, row, col) = MAT(MatA, row, col) - MAT(MatB, row, col);
    }
  }
  return matErr_None;
}

enum MatErrType MulMat(Matrix MatC, Matrix MatA, Matrix MatB,
		       int numRowsA, int numColsA,
		       int numRowsB, int numColsB) {
  
  if (!checkRowCol(MatA, numRowsA, numColsA)) {
    return matErr_Size;
  }

  if (!checkRowCol(MatB, numRowsB, numColsB)) {
    return matErr_Size;
  }

  if (!checkRowCol(MatC, numRowsA, numColsB)) {
    return matErr_Size;
  }

  if (numColsA != numRowsB) {
    return matErr_Size;
  }  
  
  for (uint32_t row = 0; row < numRowsA; row++) {
    for (uint32_t col = 0; col < numColsB; col++) {
      // Set MatC to be 0, and then perform the multiplicative sum
      MAT(MatC, row, col) = 0;
      for (uint32_t prod = 0; prod < numColsA; prod++) {
	MAT(MatC, row, col) += MAT(MatA, row, prod) * MAT(MatB, prod, col);
      }
    }
  }
  return matErr_None;
}

enum MatErrType MulMatVec(Vector VecC, Matrix MatA, Vector VecB,
			  int numRowsA, int numColsA, int numRowsB) {

  if (!checkRowCol(MatA, numRowsA, numColsA)) {
    return matErr_Size;
  }

  if (!checkIndex(VecB, numRowsB)) {
    return matErr_Size;
  }

  if (!checkIndex(VecC, numRowsA)) {
    return matErr_Size;
  }

  if (numColsA != numRowsB) {
    return matErr_Size;
  }

  for (uint32_t row = 0; row < numRowsA; row++) {
    VEC(VecC, row) = 0;
    for (uint32_t prod = 0; prod < numColsA; prod++) {
      VEC(VecC, row) += MAT(MatA, row, prod) * VEC(VecB, prod);
    }
  }
  return matErr_None;  
}

enum MatErrType ScaleMat(Matrix MatA, Matrix MatB, double scalar,
			 uint32_t numRows, uint32_t numCols) {

  if (!checkRowCol(MatA, numRows, numCols)) {
    return matErr_Size;
  }

  if (!checkRowCol(MatB, numRows, numCols)) {
    return matErr_Size;
  }
  
  if (!isValidDouble(scalar)) {
    return matErr_InvalidDouble;
  }

  for (uint32_t row = 0; row < numRows; row++) {
    for (uint32_t col = 0; col < numCols; col++) {
      MAT(B, row, col) = scalar * MAT(A, row, col);
    }
  }
  return matErr_None;
}

enum MatErrType ScaleVec(Vector VecA, Vector VecB, double scalar,
			 uint32_t numRows) {

  if (!checkIndex(VecA, numRows)) {
    return matErr_Size;
  }
  
  if (!checkIndex(VecB, numRows)) {
    return matErr_Size;
  }
  
  if (!isValidDouble(scalar)) {
    return matErr_InvalidDouble;
  }

  for (uint32_t row = 0; row < numRows; row++) {
    VEC(B, row) = scalar * VEC(A, row);
  }
  return matErr_None;
}

void TransposeMat(Matrix A, Matrix B, uint32_t numRows, uint32_t numCols) {

  for (uint32_t row = 0; row < numRows; row++) {
    for (uint32_t col = 0; col < numCols; col++) {
      MAT(B, col, row) = MAT(A, row, col);
    }
  }
}

void BackwardsSubstitution(Matrix A, Vector B, uint32_t numRows) {

  for (uint32_t i = numRows - 1; i >= 0; i--) {
    for (uint32_t j = i + 1; j < numRows; j++) {
      VEC(B, i) -= MAT(A, i, j) * VEC(B, j);
    }
    VEC(B, i) /= MAT(A, i, i);
  }
}

enum MatErrType GaussJordan(Matrix A, Vector B, uint32_t numRows) {
  
  if (!checkRowCol(A, numRows, numRows)) {
    return matErr_Size;
  }
  
  for (uint32_t k = 0; k < numRows; k++) {
    int columnMaxIndex = k; 
    /* Perform partial pivoting */
    for (int i = k + 1; i < numRows; i++) {
      double columnMaxValue = MAT(A, columnMaxIndex, k);
      double columnCandidateValue = MAT(A, i, k);
      columnMaxIndex = (fabs(columnMaxValue) < fabs(columnCandidateValue)) ?
	i : columnMaxIndex;
    }

    /* |A[pivotIndex][k]| <= 0? Absurd! */
    if (MAT(A, columnMaxIndex, k) == 0) {
      return matErr_Singular;
    }

    /* Swap row columnMaxValue */
    for (uint32_t j = 0; j < numRows; j++) {
      double* rowValueAtColumnMax = &MAT(A, columnMaxIndex, j);
      double* rowValueAtPivotIndex = &MAT(A, k, j);
      swapDouble(rowValueAtColumnMax, rowValueAtPivotIndex);
    }

    /* Swap entries in B vector */
    double* vecValAtColumnMax = &VEC(B, columnMaxIndex);
    double* vecValAtPivotIndex = &VEC(B, k);
    swapDouble(vecValAtColumnMax, vecValAtPivotIndex);

    /* Perform elimination on from row k onto all rows beneath */
    for (uint32_t i = k + 1; i < numRows; i++) {
      double eliminationScalar = MAT(A, i, k) / MAT(A, k, k);
      for (uint32_t j = k; j < numRows; j++) {
	MAT(A, i, j) -= MAT(A, k, j) * eliminationScalar; 
      }
      /* Elimination on vector B */
      VEC(B, i) -= VEC(B, k) * eliminationScalar;
    }
  }
  /* Backward substitution on B to solve for X */
  BackwardsSubstitution(A, B, numRows);
  
  return matErr_None;
}

void StandardInnerProduct(Vector A, Vector B, double* product, uint32_t numRows) {
  
  double innerProduct = 0;

  for (uint32_t k = 0; k < numRows; k++) {
    innerProduct += VEC(A, k) * VEC(B, k);
  }

  *product = innerProduct;
  
}

void l2Normalization(Vector A, double* l2norm, uint32_t numRows) {
  StandardInnerProduct(A, A, l2norm, numRows);
  *l2norm = sqrt(*l2norm);
}

enum MatErrType GramSchmidtOrthogonalization(Matrix A, Matrix Q,
					     Matrix R, int numRows,
					     uint32_t numCols) {

  Vector U, U_sub, U_proj;
  newVect(&U, numRows);
  newVect(&U_sub, numRows);
  newVect(&U_proj, numRows);

  for (uint32_t col = 0; col < numCols; col++) {
    
    /* Create U vector from col column of A */
    for (uint32_t row = 0; row < numRows; row++) {
      VEC(U, row) = MAT(A, row, col);
      VEC(U_sub, row) = VEC(U, row);
    }
    
    for (uint32_t j = 0; j < col; j++) {

      /* Get jth vector from Q which has column index less than col */
      for (uint32_t row = 0; row < numRows; row++) {
	VEC(U_proj, row) = MAT(Q, row, j);
      }

      /* Calculate projection of U onto U_proj */
      double UOntoU_proj;
      StandardInnerProduct(U, U_proj, &UOntoU_proj, numRows);
      
      /* Store projection in R */
      MAT(R, j, col) = UOntoU_proj / MAT(R, j, j);
      UOntoU_proj /= SQR(MAT(R, j, j));
      
      /* Orthogonalize */
      for (uint32_t row = 0; row < numRows; row++) {
	VEC(U_sub, row) = VEC(U_sub, row) - VEC(U_proj, row) * UOntoU_proj;
      }
    }

    double U_subl2Norm;
    l2Normalization(U_sub, &U_subl2Norm, numRows);
    MAT(R, col, col) = U_subl2Norm;
    for (uint32_t row = 0; row < numRows; row++) {
      MAT(Q, row, col) = VEC(U_sub, row);
    }    
  }

  /* Normalize Q */
  for (uint32_t col = 0; col < numCols; col++) {
    double l2NormCol = MAT(R, col, col);
    for (uint32_t row = 0; row < numRows; row++) {
      MAT(Q, row, col) = MAT(Q, row, col) / l2NormCol;
    }
  }
  
  deleteVect(&U);
  deleteVect(&U_sub);
  deleteVect(&U_proj);
  return matErr_None;
}

/* 
   QRDecomp:
     
   Parameters:
     
   Preconditions:
     
   Postconditions:
     
*/
enum MatErrType QRDecomp(Matrix A, Vector B, Matrix Q, Matrix R,
			 uint32_t numRows, uint32_t numCols) {

  if (!checkRowCol(A, numRows, numRows)) {
    return matErr_Size;
  }

  for (uint32_t i = 0; i < numCols; i++) {
    for (uint32_t j = 0; j < numCols; j++) {
      MAT(R, i, j) = 0;
    }
  }
  
  /* Perform Gram Schmidt Orthogonalization to obtain reduced Q and R. */
  enum MatErrType err =
    GramSchmidtOrthogonalization(A, Q, R, numRows, numCols);
  
  if (err != matErr_None) {
    return err;
  }

  Matrix QT;
  newMat(&QT, numCols, numRows);
  
  /* Tranpose Q and create Q^T*B */
  TransposeMat(Q, QT, numRows, numCols);

  Vector C;
  newVect(&C, numRows);

  for (uint32_t i = 0; i < numRows; i++) {
    VEC(C, i) = VEC(B, i);
  }
  
  err = MulMatVec(B, QT, C, numRows, numCols, numRows);

  if (err != matErr_None) {
    return err;
  }
  
  /* Back substitution on B using R */
  BackwardsSubstitution(R, B, numRows);
  
  deleteMat(&QT);
  deleteVect(&C);
  return matErr_None;
}

enum MatErrType GaussSeidel(Matrix A, Vector B, Vector X,
			    uint32_t numRows, uint32_t maxIter,
			    double eps1, double eps2) {

  enum opType { opContinue, opConverge,
		opSingular, opError };

  enum opType operType = opContinue;
  uint32_t iter = 0;
  Vector Xold;
  newVect(&Xold, numRows);

  /* normalize matrix A and vector B */
  for (uint32_t i = 0; i < numRows; i++) {
    double denom = MAT(A, i, i);
    if (denom < eps1) {
      return matErr_Singular;
    }
    VEC(X, i) /= denom;
    for (uint32_t j = 0; j < numRows; j++) {
      MAT(A, i, j) /= denom;
    }
  }

  /* perform Gauss-seidel iteration */
  while (operType == opContinue) {
    for (uint32_t i = 0; i < numRows; i++) {
      VEC(Xold, i) = VEC(X, i);
      VEC(X, i) = 0;
      for (uint32_t j = 0; j < numRows; j++) {
	if (j != i) {
	  VEC(X, i) -= MAT(A, i, j) * VEC(X, j);
	}
      }
      VEC(X, i) += VEC(B, i);
    }
    /* check for convergence of X */
    double devMax = fabs(VEC(Xold, 0) - VEC(X, 0)) / VEC(X, 0);
    for (uint32_t i = 1; i < numRows; i++) {
      double dev = fabs(VEC(Xold, i) - VEC(X, i)) / VEC(X, i);
      devMax = (dev > devMax) ? dev : devMax;
    }
    if (devMax <= eps2) {
      operType = opConverge;
    } else {
      iter++;
      if (iter > maxIter) {
	operType = opError;
      }
    }
  }
  deleteVect(&Xold);

  switch (operType) {
  case opConverge:
    return matErr_None;
  case opSingular:
    return matErr_Singular;
  case opError:
    return matErr_IterLimit;
  default:
    return matErr_None;
  }
}
