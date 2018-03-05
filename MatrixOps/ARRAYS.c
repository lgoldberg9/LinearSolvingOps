/* 
 * Matrix algebra package for linear system solving.
 *
 * Created by: Logan Goldberg
 * Last Modified: 3/4/2018
 *  
 */
#include "ARRAY.H"

/* Constructs a new (maxRows x maxCols) matrix */
bool newMat(Matrix* Mat, size_t maxRows, size_t maxCols) {
  maxRows = (maxRows < 1) ? 1 : maxRows;
  maxCols = (maxCols < 1) ? 1 : maxCols;

  Mat->maxRows = maxRows;
  Mat->maxCols = maxCols;
  Mat->pData = (double*) malloc(maxRows * maxCols * sizeof(double));
  
  return (Mat->pData != NULL);
}

/* Frees Mat */
void deleteMat(Matrix* Mat) {
  if (Mat->pData) {
    free(Mat->pData);
  }
  Mat->maxCols = 0;
  Mat->maxRows = 0;
}

/* Generates new double vector on parameter Vect */
bool newVect(Vector* Vect, size_t maxSize) {
  maxSize = (maxSize < 1) ? 1 : maxSize;
  Vect->maxSize = maxSize;
  Vect->pData = (double*) malloc(maxSize * sizeof(double));
  return (Vect->pData != NULL);
}

/* Frees Vect */
void deleteVect(Vector* Vect) {
  if (Vect->pData) {
    free(Vect->pData);
  }
  Vect->maxSize = 0;
}

/* Checks if row and col are within bounds */
bool checkRowCol(Matrix Mat, uint32_t row, uint32_t col) {
  return (row > 0 && col > 0 &&
	  row <= Mat.maxRows && col <= Mat.maxCols);
}

/* Check if maxSize is within vect bounds */
bool checkIndex(Vector Vect, size_t maxSize) {
  return (maxSize > 0 && maxSize <= Vect.maxSize);
}

bool isValidDouble(double x) {
  return (x != INFINITY && x != NAN);
}
