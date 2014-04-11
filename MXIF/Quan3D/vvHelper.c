/** \file vvHelper.c
\brief This is the implementation file for helper functions used in CCPI 
VolView plug ins

\author David Worth STFC
\date May 2012
*/

#include <stdlib.h>

double **CovMat(double **matrix, int numRows, int numCols) {

  /* Allocate memory for the returned matrix */
  double* *covMatrix = malloc(sizeof(double*)*numCols);
  /* Mean values */
  double *meanVals = malloc(sizeof(double)*numCols);
  int i, j, k;
  /* Pointer to rows of input matrix */
  double *row = NULL;
  
  /* Find the mean values for the columns */
  for (j = 0; j < numCols; j++) meanVals[j] = 0.0;
  for (i = 0; i < numRows; i++) {
    row = matrix[i];
    for (j = 0; j < numCols; j++) {
      meanVals[j] += row[j];
    }
  }
  for (j = 0; j < numCols; j++) meanVals[j] /= numRows;
  
  /* Calculate the covariance matrix */
  for (j = 0; j < numCols; j++) {
    covMatrix[j] = malloc(sizeof(double)*numCols);
    for (i = 0; i < numCols; i++) covMatrix[j][i] = 0.0;
  }
  
  for (i = 0; i < numRows; i++) {
    for (j = 0; j < numCols; j++) {
      for (k = 0; k < numCols; k++) {
        covMatrix[j][k] += ((matrix[i][j]-meanVals[j]) * 
                            (matrix[i][k]-meanVals[k]))/(double)(numRows-1);
      }
    }
  }
  
  return covMatrix;
}
