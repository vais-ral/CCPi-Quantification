/** \file vvHelper.h
\brief This is the header file for helper functions used in CCPI VolView plug
ins

\author David Worth STFC
\date May 2012
*/

/**
\brief Calculate the covariance matrix for the given matrix

This function will allocate the memory for the returned matrix <b>it is up to
the calling function to free this memory.</b>

\param matrix The matrix to find the covariance matrix for
\param numRows The number of rows in the matrix
\param numCols The number of columns in the matrix
\return The numCols x numCols covariance matrix
*/
double **CovMat(double **matrix, int numRows, int numCols);
