/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file Quan3D.cpp
This is the implementation file for class CCPiQuantification3D<IT>
*/

#include "Quan3D.hpp"
#include "QuanWorker.hpp"

#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "itkArray.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkVariableSizeMatrix.h"
#include "vtkCellArray.h"
#include "vtkImageImport.h"
#include "vtkMassProperties.h"
#include "vtkMath.h"
#include "vtkParametricEllipsoid.h"
#include "vtkPointsProjectedHull.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkSynchronizedTemplates3D.h"

/*
 * Destructor. Cleans up allocated memory.
 */
template <class IT>
CCPiQuantification3D<IT>::~CCPiQuantification3D()
{
  if (m_WrappedImage != NULL) delete[] m_WrappedImage;
  if (m_ImageOrigin != NULL) delete[] m_ImageOrigin;
  if (m_ImageSpacing != NULL) delete[] m_ImageSpacing;
}

/*
 * Initialise class with separate values
 *
 * \param imageData Pointer to raw image data
 * \param dimensions Number of voxels in each direction (x,y,z)
 * \param numComponents Number of values for each voxel
 * \param inputScalarType Type of data in the image
 * \param origin The position of the origin for this data (x,y,z)
 * \param spacing The physical spacing of the voxels in each direction
 * \param minValue The minimum value of the data
 * \param maxValue The maximum value of the data
 */
template <class IT>
void CCPiQuantification3D<IT>::Initialise(IT *imageData, const int *dimensions, 
  int numComponents, int inputScalarType, const float *origin, const float *spacing, 
  float minValue, float maxValue)
{
  m_ImageData = imageData;
  m_Dimensions = (int*)dimensions;
  m_NumComponents = numComponents;
  m_InputVolumeScalarType = inputScalarType;
  m_ImageOrigin = new float[3];
  m_ImageOrigin[0] = origin[0];
  m_ImageOrigin[1] = origin[1];
  m_ImageOrigin[2] = origin[2];
  m_ImageSpacing = new float[3];
  m_ImageSpacing[0] = spacing[0];
  m_ImageSpacing[1] = spacing[1];
  m_ImageSpacing[2] = spacing[2];
  m_MinValue = minValue, 
  m_MaxValue = maxValue;
}

/*
 * Set the physical voxel size so output has physical size
 * \param voxelSize Physical voxel size
 */
template <class IT>
void CCPiQuantification3D<IT>::SetVoxelSize(float voxelSize)
{
  m_VoxelSize = voxelSize;
  m_AreaFactor = pow(voxelSize,2);
  m_VolumeFactor = pow(voxelSize,3);
}

/*
 * Create the list of voxel values and list of x,y,z indices for each value.
 * This is a necessary set up step before any real calculation can be 
 * performed
 */
template <class IT>
void CCPiQuantification3D<IT>::CreateVoxelIndexList(void)
{
  /* 
   * Create array of voxel counts for each voxel value except the min value 
   * which we take to be the background. Also count the number of different
   * voxel values and all the non-background voxels.
   */
  m_MaxVoxelVals = (int)ceil(m_MaxValue)-(int)floor(m_MinValue)+1;
  m_VoxelCounts.SetSize(m_MaxVoxelVals);
  m_VoxelCounts.Fill(0);

  m_TotalParticles = 0; 
  m_SumVoxels = 0;

  for (int i = 0; i < m_Dimensions[0]*m_Dimensions[1]*m_Dimensions[2]; i++) {
    if (m_ImageData[i] != m_MinValue) { /* Not the minimum value */
      /* Not seen this value before so it's a new particle */
      if (m_VoxelCounts[m_ImageData[i]] == 0) m_TotalParticles++;
      m_VoxelCounts[m_ImageData[i]]++;
      m_SumVoxels++;
    }
  }

  /* 
   * Create arrays that will hold the (i,j,k) coordinates of each voxel with
   * a given intensity value. Each of these arrays will be stored in another
   * array - the index of which is the intensity value for those voxels
   */
  m_VoxelArray.reserve(m_MaxVoxelVals);
  std::vector<int> indexVoxelSubArray(m_MaxVoxelVals,0);

  for (int i = 0; i < m_MaxVoxelVals; i++) {
    std::vector< itk::FixedArray<int,3> > voxelSubArray(m_VoxelCounts[i]);
    m_VoxelArray.push_back(voxelSubArray);
  }

  /* Now populate that array with values */
  std::vector< itk::FixedArray<int,3> > voxelSubArray;
  itk::FixedArray<int,3> coordinates;
  /* Don't change m_ImageData so use a separate pointer */
  IT *inPtr = m_ImageData; 
  /* loop over the slices */
  for (int k = 0; k < m_Dimensions[2]; k++ ) {                                      
    /* loop over the rows and handle aborts */
    for (int j = 0; j < m_Dimensions[1]; j++ ) {                         
      /* loop over the columns */
      for (int i = 0; i < m_Dimensions[0]; i++ ) {
        if (*inPtr != m_MinValue) {
          /* This is a bit complicated but saves creating temporary objects */
          ((m_VoxelArray[*inPtr])[ indexVoxelSubArray[*inPtr] ])[0] = i;
          ((m_VoxelArray[*inPtr])[ indexVoxelSubArray[*inPtr] ])[1] = j;
          ((m_VoxelArray[*inPtr])[ indexVoxelSubArray[*inPtr] ])[2] = k;
          indexVoxelSubArray[*inPtr]++;
        }
        inPtr++;
      }
    }
  }         
}

/*
 * Prepare to iterate through the voxel values. Must be called before
 * NextValueQuantification().
 */
template <class IT>
void CCPiQuantification3D<IT>::PrepareForQuantification()
{
  // Create and intialise image to be used in iso-surface construction.
  // Could be big so only do it once.
  m_WrappedImage = new IT[(m_Dimensions[0]+2)*(m_Dimensions[1]+2)*
    (m_Dimensions[2]+2)];
  memset(m_WrappedImage, (IT)0, (m_Dimensions[0]+2)*(m_Dimensions[1]+2)*
    (m_Dimensions[2]+2)*sizeof(IT));
/*  m_WrappedImage = new IT[m_Dimensions[0]*m_Dimensions[1]*m_Dimensions[2]];
  memset(m_WrappedImage, (IT)0, m_Dimensions[0]*m_Dimensions[1]*
    m_Dimensions[2]*sizeof(IT));*/

  // Initialise counter of next voxel value to -1 since GetNextVoxelValue must
  // increment value before it starts to avoid getting stuck
  m_NextValue = -1;

  // Write header for CSV values
  WriteCSVHeader();

}

/**
 * Get the object to calculate quantification data for object index i
 *
 * \param i The index of the calculation object required
 * \return Pointer to calculation object which the caller should destroy.
 */
template <class IT>
CCPiQuantificationWorker<IT>* CCPiQuantification3D<IT>::GetWorker(int i)
{
  return new CCPiQuantificationWorker<IT>(this, i);
}

/**
 * Get the next object to calculate quantification data 
 *
 * \return Pointer to calculation object which the caller should destroy.
 */
template <class IT>
CCPiQuantificationWorker<IT>* CCPiQuantification3D<IT>::GetNextWorker()
{
  int value = GetNextVoxelValue();
  if (value ==  m_MaxVoxelVals) {
    return NULL;
  }
  else {
    return new CCPiQuantificationWorker<IT>(this, value);
  }
}



/*
 * Print summary data of the volume to cout
 */
template <class IT>
void CCPiQuantification3D<IT>::PrintSummaryData(void)
{
  /* Print out some data about the image */
  cout << "The dimensions of the image are: " << m_Dimensions[0] << " " << m_Dimensions[1] << 
          " " << m_Dimensions[2] << endl;
  cout << "The scalar data has size " << sizeof(IT) << endl;

  cout << "There are " << m_NumComponents << " components" << endl;
  cout << "The minimum value in the image is " << m_MinValue;
  cout << " and the maximum value is " << m_MaxValue << endl;

  double porosity = (double)m_SumVoxels/(double)(m_Dimensions[0]*
    m_Dimensions[1]*m_Dimensions[2]);
  cout << "Porosity is " << porosity << endl;
  cout << "Number of particles is " << m_TotalParticles << endl;
}

/**
 * Write headers for comma separated variable data table to data stream
 */
template <class IT>
void CCPiQuantification3D<IT>::WriteCSVHeader(void)
{
  m_FinalStatisticsLog << "Label, Border, Number_of_Voxels, minX, minY, minZ, maxX, maxY, maxZ, ";
  m_FinalStatisticsLog << "meanX, meanY, meanZ, Bbox_diag, VoxelCountsVolume, VoxelCountDia, ";
  m_FinalStatisticsLog << "PCA_2D_area, PCA_2D_dia, Surface_Area, Surface_Volume, Surface_Dia, ";
  m_FinalStatisticsLog << "Sphericity, RSaVol" << endl; 
}

/**
 * Write statistics calculated so far to file.
 *
 * \param filename Name of file to print to.
 */
template <class IT>
void CCPiQuantification3D<IT>::WriteCSVData(std::string filename)
{
  // Create a file and write to it
  std::ofstream finalStatisticsLog(filename.c_str(), std::ios::trunc);
  finalStatisticsLog << m_FinalStatisticsLog.str().c_str();
  finalStatisticsLog.close();
}


/*
\brief Calculate sphere diameter from the sphere volume.

Use the formula \f$ 2*(3v/(4\pi))^{\frac{1}{3}} \f$ where \f$ v \f$ is the 
sphere volume.

\param volume The sphere volume
\return Diameter of the sphere
*/
template <class IT>
float CCPiQuantification3D<IT>::Vol2Dia(float Volume)
{
#if VTK_MAJOR_VERSION > 5 
  return 2.0 * pow((3.0*Volume)/(4.0*vtkMath::Pi()),(1.0/3.0));
#else
  return 2.0 * pow((3.0*Volume)/(4.0*vtkMath::DoublePi()),(1.0/3.0));
#endif
}

/*
\brief Calculate the maximum, minimum and mean positions from a set of (i,j,k)
coordinates.

Use itk::Vector for return objects to make subsequent calculations with them
more straightforward.

\param Positions Reference to positions we are working with
\param Minimum Minimum values for (i,j,k) from positions
\param Maximum Maximum values for (i,j,k) from positions
\param Mean Mean values for (i,j,k) from positions
*/
template <class IT>
void CCPiQuantification3D<IT>::GetPositionStats(
  const std::vector< itk::FixedArray<int,3> >& Positions, 
  itk::Vector<int,3>& Minimum, itk::Vector<int,3>& Maximum, 
  itk::Vector<double,3>& Mean)
{
  /* Initialise return data */
  Minimum.Fill(INT_MAX);
  Maximum.Fill(INT_MIN);
  Mean.Fill(0);

  /* Iterate over positions and find min, max, mean of the individual 
   * components */
  std::vector< itk::FixedArray<int,3> >::const_iterator iter;
  for (iter = Positions.begin(); iter < Positions.end(); iter++) {
    if ( (*iter)[0] < Minimum[0] ) Minimum[0] = (*iter)[0];
    if ( (*iter)[0] > Maximum[0] ) Maximum[0] = (*iter)[0];
    Mean[0] += (*iter)[0];
    if ( (*iter)[1] < Minimum[1] ) Minimum[1] = (*iter)[1];
    if ( (*iter)[1] > Maximum[1] ) Maximum[1] = (*iter)[1];
    Mean[1] += (*iter)[1];
    if ( (*iter)[2] < Minimum[2] ) Minimum[2] = (*iter)[2];
    if ( (*iter)[2] > Maximum[2] ) Maximum[2] = (*iter)[2];
    Mean[2] += (*iter)[2];
  }

  Mean /= Positions.size();
}

/*
\brief Calculate the covariance matrix from the given array of (i,j,k) 
coordinates

\param positions Reference to positions we are working with
\param covarianceMatrix The calculated covariance matrix
*/
template <class IT>
void CCPiQuantification3D<IT>::CalculateCovarianceMatrix(
  const std::vector< itk::FixedArray<int,3> >& Positions, 
  itk::Matrix<double,3,3>& CovarianceMatrix)
{
  itk::Vector<double,3> mean;
  mean.Fill(0);
  
  /* Iterate over positions and find min, max, mean of the individual 
   * components */
  std::vector< itk::FixedArray<int,3> >::const_iterator iter;
  for (iter = Positions.begin(); iter < Positions.end(); iter++) {
    mean[0] += (*iter)[0];
    mean[1] += (*iter)[1];
    mean[2] += (*iter)[2];
  }
  mean /= Positions.size();

  /* Calculate the covariance matrix */
  CovarianceMatrix.Fill(0);

  for (iter = Positions.begin(); iter < Positions.end(); iter++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        CovarianceMatrix[j][k] += (((*iter)[j]-mean[j]) *
                            ((*iter)[k]-mean[k]))/(double)(Positions.size()-1);
      }
    }
  }
/*  std::cout << "Covariance matrix\n";
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      std::cout << CovarianceMatrix[j][k] << " ";
    }
    std::cout << std::endl;
  }*/
}

/**
 * Work out which is the next voxel value to work on.
 * Updates m_NextValue and returns that value for good measure!
 *
 * \return The voxel value to work on next
 */
template <class IT>
int CCPiQuantification3D<IT>::GetNextVoxelValue() {
  
  // Incremenet the value to avoid getting stuck
  m_NextValue++;
  while ((m_NextValue < m_MaxVoxelVals) && (m_VoxelCounts[m_NextValue] == 0)) 
    m_NextValue++;

  return m_NextValue;
}

template <class IT>
void CCPiQuantification3D<IT>::SetQuantificationResultByWorker(int labelIndex, std::map<std::string, double> result)
{
	for(std::map<std::string,double>::iterator itr=result.begin(); itr!= result.end();itr++)
		QuantificationResult.SetValue(itr->first,labelIndex,itr->second);
}

template class CCPiQuantification3D<double>;
template class CCPiQuantification3D<float>;
template class CCPiQuantification3D<long>;
template class CCPiQuantification3D<long long>;
template class CCPiQuantification3D<unsigned long>;
template class CCPiQuantification3D<unsigned long long>;
template class CCPiQuantification3D<short>;
template class CCPiQuantification3D<unsigned short>;
template class CCPiQuantification3D<char>;
template class CCPiQuantification3D<unsigned char>;
template class CCPiQuantification3D<signed char>;
template class CCPiQuantification3D<int>;
template class CCPiQuantification3D<unsigned int>;

