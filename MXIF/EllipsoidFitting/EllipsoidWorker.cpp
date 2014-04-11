/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file EllipsoidWorker.cpp
This is the implementation file for class CCPiEllipsoidWorker<IT>
*/

#include "EllipsoidWorker.hpp"

#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "itkArray.h"
#include "itkEllipsoidInteriorExteriorSpatialFunction.h"
#include "itkPoint.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkVariableSizeMatrix.h"
#include "vtkImageData.h"
#include "vtkImageImport.h"


/**
 * Destructor. Cleans up allocated memory.
 */
template <class IT>
CCPiEllipsoidWorker<IT>::~CCPiEllipsoidWorker()
{
  m_FittingData = NULL;
}


/**
 * Do the calculation
 * \return 0 if calculation was run or 1 if this value had insufficient values
 * to be significant.
 */
template <class IT> 
int CCPiEllipsoidWorker<IT>::Run()
{
  int voxelCount;
  itk::Vector<int,3> maximum, minimum;
  itk::Vector<double,3> mean;
  itk::Matrix<double,3,3> covarianceMatrix;
  itk::SymmetricEigenAnalysis< itk::Matrix<double,3,3>, itk::Vector<double,3>,
                               itk::Matrix<double,3,3> > eigenSystem(3);
  itk::Vector<double,3> eigenValues;
  itk::Matrix<double,3,3> eigenVectors;

  // Get data for this region from the image data
  voxelCount = m_FittingData->m_VoxelCounts[m_Id];
  // If this feature doesn't have enough voxels to be significant do nothing
  if ((voxelCount == 0) || (voxelCount < m_FittingData->m_MinFeatureSize)) {
    return 1;
  }
  std::vector< itk::FixedArray<int,3> > voxelPositions = m_FittingData->m_VoxelArray[m_Id];
  m_FittingData->GetPositionStats(voxelPositions, minimum, maximum, mean);

  // Calculate the covariance matrix
  m_FittingData->CalculateCovarianceMatrix(voxelPositions, covarianceMatrix);

  // Calculate the eigenvalues and vectors
  eigenSystem.ComputeEigenValuesAndVectors(covarianceMatrix, eigenValues,
                                           eigenVectors);

  // Store ellipsoid parameters
  for (int i = 0; i < 3; i++) {
    m_SemiRadii[i] = sqrt(eigenValues[i])*2.0;
    m_Centre[i] = mean[i];
    for (int j = 0; j < 3; j++) m_Orientation[i*3+j] = eigenVectors(i,j);
  }
  
  // Write first lot of data for output file
  m_FinalStatisticsLog << m_Id << ", " << sqrt(eigenValues[0])*2.0 << ", ";
  m_FinalStatisticsLog << eigenVectors(0,0) << ", " << eigenVectors(0,1) << ", " << eigenVectors(0,2) << ", ";
  m_FinalStatisticsLog << sqrt(eigenValues[1])*2.0 << ", ";
  m_FinalStatisticsLog << eigenVectors(1,0) << ", " << eigenVectors(1,1) << ", " << eigenVectors(1,2) << ", ";
  m_FinalStatisticsLog << sqrt(eigenValues[2])*2.0 << ", ";
  m_FinalStatisticsLog << eigenVectors(2,0) << ", " << eigenVectors(2,1) << ", " << eigenVectors(2,2) << ", ";

  // Write bounding box to output file
  m_FinalStatisticsLog << minimum[0] << ", " << minimum[1] << ", " << minimum[2] << ", ";
  m_FinalStatisticsLog << maximum[0] << ", " << maximum[1] << ", " << maximum[2] << ", ";
  m_FinalStatisticsLog << mean[0] << ", " << mean[1] << ", " << mean[2] << std::endl;


  return 0;

}


/**
 * Draw ellipsoid for this data
 * 
 * \param outputImage Output image raw data array
 */
template <class IT> 
void CCPiEllipsoidWorker<IT>::DrawEllipsoid(unsigned short *outputImage)
{

  // If there is somewhere to output the ellipsoid to the do so
  if (outputImage != NULL) {
    itk::EllipsoidInteriorExteriorSpatialFunction<3>::Pointer ellipsoid = 
      itk::EllipsoidInteriorExteriorSpatialFunction<3>::New();
      
    // Set the axes lengths
    itk::EllipsoidInteriorExteriorSpatialFunction<3>::InputType axes;
    axes[0] = 2.0*m_SemiRadii[0];
    axes[1] = 2.0*m_SemiRadii[1];
    axes[2] = 2.0*m_SemiRadii[2];
    //std::cout << axes[0] << ", " << axes[1] << ", " << axes[2] << std::endl;
    ellipsoid->SetAxes(axes);
     
    // Set the centre
    itk::EllipsoidInteriorExteriorSpatialFunction<3>::InputType centre;
    centre[0] = m_Centre[0];
    centre[1] = m_Centre[1];
    centre[2] = m_Centre[2];
    ellipsoid->SetCenter(centre);
    //std::cout << m_Centre[0] << ", " << m_Centre[1] << ", " << m_Centre[2] << std::endl;
   
    // Set the axes orientations
    vnl_matrix<double> orientations(m_Orientation, 3, 3);
    ellipsoid->SetOrientations(orientations);
/*    std::cout << m_Orientation[0] << ", " << m_Orientation[1] << ", " << m_Orientation[2] << std::endl;
    std::cout << m_Orientation[3] << ", " << m_Orientation[4] << ", " << m_Orientation[5] << std::endl;
    std::cout << m_Orientation[6] << ", " << m_Orientation[7] << ", " << m_Orientation[8] << std::endl;*/
   
    // Now iterate through the output image and set those voxels in the 
    // ellipsoid to the id of this set of voxels and leave others alone
    double testPosition[3];
    int index = 0;
    for (int i = 0; i < m_FittingData->m_Dimensions[0]; i++) {
      for (int j = 0; j < m_FittingData->m_Dimensions[1]; j++) {
        for (int k = 0; k < m_FittingData->m_Dimensions[2]; k++) {
          testPosition[0] = i;
          testPosition[1] = j;
          testPosition[2] = k;
          index = k*m_FittingData->m_Dimensions[0]*m_FittingData->m_Dimensions[1] + 
              j*m_FittingData->m_Dimensions[0]+i;
          if (ellipsoid->Evaluate(testPosition) == 1) {
            outputImage[index] = (unsigned short)m_Id;
          }
        }
      }
    }

	  //ellipsoid->Delete();
  }
}
 


/**
 * Write statistics calculated so far to file.
 *
 * \param filename Name of file to print to.
 */
template <class IT>
void CCPiEllipsoidWorker<IT>::WriteCSVData(std::string filename)
{
  // Create a file and write to it
  std::ofstream finalStatisticsLog(filename.c_str(), std::ios::app);
  finalStatisticsLog << m_FinalStatisticsLog.str();
  finalStatisticsLog.close();
}


template class CCPiEllipsoidWorker<double>;
template class CCPiEllipsoidWorker<float>;
template class CCPiEllipsoidWorker<long>;
template class CCPiEllipsoidWorker<long long>;
template class CCPiEllipsoidWorker<unsigned long>;
template class CCPiEllipsoidWorker<unsigned long long>;
template class CCPiEllipsoidWorker<short>;
template class CCPiEllipsoidWorker<unsigned short>;
template class CCPiEllipsoidWorker<char>;
template class CCPiEllipsoidWorker<unsigned char>;
template class CCPiEllipsoidWorker<signed char>;
template class CCPiEllipsoidWorker<int>;
template class CCPiEllipsoidWorker<unsigned int>;

