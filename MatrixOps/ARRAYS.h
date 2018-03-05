/* 
 * Matrix algebra package for linear system solving.
 *
 * Created by: Logan Goldberg
 * Last Modified: 3/4/2018
 *  
 */

#ifndef __ARRAYS_H__
#define __ARRAYS_H__

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#define MAT(Mat, row, col) (Mat.pData[(row) + Mat.maxRows * (col)])
#define VEC(Vect, index) (Vect.pData[(index)])

typedef struct VectTag {
  double* pData;
  size_t maxSize;
} Vector;

typedef struct MatTag {
  double* pData;
  size_t maxRows;
  size_t maxCols;
} Matrix;

/* 
   Indexing is as follows: if m is the ijth entry of matrix C,
   then m = c->pData[i + c.maxRows * j] is the jith entry. 
*/

/* Generates new matrix on parameter Mat
 * Returns whether matrix initialized
 */
bool newMat(Matrix* Mat, size_t maxRows, size_t maxCols);

/* Frees Mat */
void deleteMat(Matrix* Mat);

/* Generates new double vector on parameter Vect 
 * Returns whether vector initialized
 */
bool newVect(Vector* Vect, size_t maxSize);

/* Frees Vect */
void deleteVect(Vector* Vect);

/* Checks if row and col are within bounds */
bool checkRowCol(Matrix Mat, uint32_t row, uint32_t col);

/* Check if maxSize is within vect bounds */
bool checkIndex(Vector Vect, size_t maxSize);

/* Checks if double is NaN or infinity */
bool isValidDouble(double x);

#endif
