/* 
   Copyright (c) 1995 Namir C. Shammas
   Version 1.0, Date 8/9/1994
   C module for basic vector and array operations:
   + Add matrices
   + Subtract matrices
   + Multiply matrices
   + Solve a set of linear equations using the
           Gauss-Jordan method
   + Solve a set of linear equations using the 
           QR decomposition method
   + Solve a set of linear equations using the
           Gauss-Seidel method
   Last modified by Logan Goldberg on 3/4/2018
*/

#ifndef __MATRIXOPS_H__
#define __MATRIXOPS_H__

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "GLOBAL.h"
#include "ARRAYS.h"

#define MATRIXOPS_EPSILON 1.0e-15
#define MATRIXOPS_BAD_RESULT -1.0e+30

enum MatErrType { matErr_None, matErr_Size, matErr_Singular,
		  matErr_IllConditioned, matErr_IterLimit,
		  matErr_InvalidDouble };

/* C module for basic vector and array operations. */
void swapInt(uint32_t* i1, uint32_t* i2);

void swapDouble(double* d1, double* d2);

enum MatErrType CopyMat(Matrix MatB, Matrix MatA,
			uint32_t numRows, uint32_t numCols);

enum MatErrType AddMat(Matrix MatC, Matrix MatA, Matrix MatB,
		       uint32_t numRows, uint32_t numCols);

enum MatErrType SubMat(Matrix MatC, Matrix MatA, Matrix MatB,
		       uint32_t numRows, uint32_t numCols);

enum MatErrType MulMat(Matrix MatC, Matrix MatA, Matrix MatB,
		       uint32_t numRowsA, uint32_t numColsA,
		       uint32_t numRowsB, uint32_t numColsB);

enum MatErrType ScaleMat(Matrix MatA, Matrix MatB, double scalar,
			 uint32_t numRows, uint32_t numCols);

enum MatErrType ScaleVec(Vector VecA, Vector VecB, double scalar,
			 uint32_t numRows);

/* 
   MulMatVec:
     Because types in C are annoying. Performs the same operation as MulMat
     except that the inputs are two vectors and a matrix.
   Parameters:
     VecC, a vector
     VecB, a vector
     MatA, a matrix
     numRowsA, number of rows in matrix A
     numColsA, number of columns in matrix A
     numRowsB, number of rows in vector B
   Preconditions:
     VecC, MatA, VecB must all be allocated.
     numRowsA, numColsA, numRowsB > 0
   Postconditions:
     VecC is the product of matrix A acting on vector B as
     AB = C.
 */
enum MatErrType MulMatVec(Vector VecC, Matrix MatA, Vector VecB,
			  uint32_t numRowsA, uint32_t numColsA,
			  uint32_t numRowsB);

/* 
   MatTranspose:
     Tranposes a matrix.
   Parameters:
     A, a matrix
     B, a matrix
     numRows, rows of matrix A
     numCols, columns of matrix A
   Preconditions:
     A must be malloced.
     numRows, numCols > 0
     A.numRows = B.numCols
     B.numRows = A.numCols
   Postconditions:
     B = A^T or in other words, B is the transpose of A and has
     dimensions numCols x numRows if A is numRows by numCols.
*/
void TransposeMat(Matrix A, Matrix B, uint32_t numRows, uint32_t numCols);

/* 
   BackwardsSubstitution:
     Solves an upper triangular matrix for a vector X given a vector B.
   
   Parameters:
     A, an upper triangular matrix
     B, a vector
     numRows, an integer
   Preconditions:
     numRows of A = numRows of B
     numRows > 0
   Postconditions:
     B is the vector of coefficients such that A multiplied from the left is B.
*/
void BackwardsSubstitution(Matrix A, Vector B, uint32_t numRows);

/* 
   GaussJordan:
     Performs Gauss-Jordan elimination on a matrix A.
     Then uses backward substitution to solve for the vector X
     corresponding to the vector B.
   Parameters:
     A, a matrix
     B, a vector
     numRows, an integer signifying the size of A.
   Preconditions:
     A is nonsingular
     dim(row(A)) = dim(row(B))
   Postconditions:
     If successful, return matErr_None
     Otherwise, return an error, most likely because A is singular.
*/
enum MatErrType GaussJordan(Matrix A, Vector B, int numRows);

/* 
   StandardInnerProduct:
     Compute the standard inner product of vector A with vector B.
   Parameters:
     A, a vector
     B, a vector
     proj, a double to store the inner product result within
     numRows, the number of entries in each vector
   Preconditions:
     A and B should be malloced.
     proj is not a null pointer.
     numRows > 0
   Postconditions:
     proj* = sum_(i=0)^(numRows - 1) A[i]*B[i]
   
*/
void StandardInnerProduct(Vector A, Vector B, double* product, uint32_t numRows);

/* 
   l2Normalization:
     Normalize a vector by its l2 norm.
   Parameters:
     A, a vector
     l2norm, the l2norm of A
     numRows, number of rows in A
   Preconditions:
     A must be malloced
     numRows > 0
   Postconditions:
     l2norm is the l2 norm of A
*/
void l2Normalization(Vector A, double* l2norm, uint32_t numRows);

/*
  
 */
enum MatErrType GramSchmidtOrthogonalization(Matrix A, Matrix Q,
					     Matrix R, uint32_t numRows,
					     uint32_t numCols);

/* 
   QRDecomp:
     
   Parameters:
     
   Preconditions:
   
   Postconditions:
     
*/
enum MatErrType QRDecomp(Matrix A, Vector B, Matrix Q, Matrix R,
			 uint32_t numRows, uint32_t numCols);

enum MatErrType GaussSeidel(Matrix A, Vector B, Vector X,
			    uint32_t numRows, uint32_t maxIter,
			    double eps1, double eps2);
#endif
