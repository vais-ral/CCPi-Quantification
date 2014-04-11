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

#include "vtkVVPluginAPI.h"

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
 * Initialise the class with the plugin information and image data passed
 * to a VolView plugin.
 * 
 * \param info Pointer to a vtkVVPluginInfo object
 * \param pds Pointer to a vtkVVProcessDataStruct object
 */
template <class IT>
void CCPiQuantification3D<IT>::Initialise(void *info, void *pds)
{
  m_ImageData = (IT *)((vtkVVProcessDataStruct*)pds)->inData;
  m_Dimensions = ((vtkVVPluginInfo*)info)->InputVolumeDimensions;
  m_NumComponents = ((vtkVVPluginInfo*)info)->InputVolumeNumberOfComponents;
  m_InputVolumeScalarType = ((vtkVVPluginInfo*)info)->InputVolumeScalarType;
  // Use double for image origin and spacing in this class but values from info
  // are float so have to do some fiddling!
  m_ImageOrigin = new double[3];
  m_ImageOrigin[0] = ((vtkVVPluginInfo*)info)->InputVolumeOrigin[0];
  m_ImageOrigin[1] = ((vtkVVPluginInfo*)info)->InputVolumeOrigin[1];
  m_ImageOrigin[2] = ((vtkVVPluginInfo*)info)->InputVolumeOrigin[2];
  m_ImageSpacing = new double[3];
  m_ImageSpacing[0] = ((vtkVVPluginInfo*)info)->InputVolumeSpacing[0];
  m_ImageSpacing[1] = ((vtkVVPluginInfo*)info)->InputVolumeSpacing[1];
  m_ImageSpacing[2] = ((vtkVVPluginInfo*)info)->InputVolumeSpacing[2];
  m_MinValue = ((vtkVVPluginInfo*)info)->InputVolumeScalarRange[0]; 
  m_MaxValue = ((vtkVVPluginInfo*)info)->InputVolumeScalarRange[1];
  // Voxel size is taken from image spacing in x-direction
  SetVoxelSize(((vtkVVPluginInfo*)info)->InputVolumeSpacing[0]);
  m_MinFeatureSize = atof(((vtkVVPluginInfo*)info)->GetGUIProperty(info, 0, VVP_GUI_VALUE));
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
void CCPiQuantification3D<IT>::Initialise(IT *imageData, int *dimensions, 
  int numComponents, int inputScalarType, double *origin, double *spacing, 
  double minValue, double maxValue)
{
  m_ImageData = imageData;
  m_Dimensions = dimensions;
  m_NumComponents = numComponents;
  m_InputVolumeScalarType = inputScalarType;
  m_ImageOrigin = new double[3];
  m_ImageOrigin[0] = origin[0];
  m_ImageOrigin[1] = origin[1];
  m_ImageOrigin[2] = origin[2];
  m_ImageSpacing = new double[3];
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
void CCPiQuantification3D<IT>::SetVoxelSize(double voxelSize)
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
 * Do quantification for the next voxel value (internal variable).
 *
 * \return 0 when there are no more values to do
 */
template <class IT>
int CCPiQuantification3D<IT>::NextValueQuantification()
{
  int voxelCount;
  double voxelCountVolume, voxelCountDia;
  itk::Vector<int,3> maximum, minimum, boundingBoxSize;
  itk::Vector<double,3> mean;
  itk::Matrix<double,3,3> covarianceMatrix;
  itk::SymmetricEigenAnalysis< itk::Matrix<double,3,3>, itk::Vector<double,3>,
                               itk::Matrix<double,3,3> > eigenSystem(3);
  itk::Vector<double,3> eigenValues;
  itk::Matrix<double,3,3> eigenVectors;

  itk::VariableSizeMatrix<double> positions;
  int numAzimuthal = 12; /* Azimuthal angle from 0 to 2PI */
  int numPolar = 6; /* Polar angle from 0 to PI */

  itk::VariableSizeMatrix<double> ellipsoidVertices;

  /* Check if we've finished */
  if (GetNextVoxelValue() ==  m_MaxVoxelVals) return 0;

  /* Do some simple stuff */
  voxelCount = m_VoxelCounts[m_NextValue];
  // If this feature doesn't have enough voxels to be significant do nothing
  if ((voxelCount == 0) || (voxelCount < m_MinFeatureSize)) {
    m_NextValue++;
    return 1;
  }

  voxelCountVolume = voxelCount * m_VolumeFactor;
  voxelCountDia = m_VoxelSize*Vol2Dia(voxelCount);

  // Start writing CSV data
  m_FinalStatisticsLog << m_NextValue << ", ";

  std::vector< itk::FixedArray<int,3> > voxelPositions = m_VoxelArray[m_NextValue];
  GetPositionStats(voxelPositions, minimum, maximum, mean);

  boundingBoxSize = maximum - minimum + 1;
  double boundingBoxDiam = boundingBoxSize.GetNorm()*m_VoxelSize;
  
  /* Calculate the covariance matrix */  
  CalculateCovarianceMatrix(voxelPositions, covarianceMatrix);

  /* Calculate the eigenvalues and vectors */
  eigenSystem.ComputeEigenValuesAndVectors(covarianceMatrix, eigenValues,
                                           eigenVectors);

  /* Calculate some new positions */
  positions.SetSize(voxelCount,3);
  int numBoundaryVoxels = 0;
  for (int k = 0; k < voxelCount; k++) {
      positions[k][0] = voxelPositions[k][0];
      if ( (voxelPositions[k][0] == 0) || (voxelPositions[k][0] == m_Dimensions[0]-1) ||
           (voxelPositions[k][1] == 0) || (voxelPositions[k][1] == m_Dimensions[2]-1) ||
           (voxelPositions[k][2] == 0) || (voxelPositions[k][2] == m_Dimensions[2]-1) )
        numBoundaryVoxels++;
      positions[k][1] = voxelPositions[k][1];
      positions[k][2] = voxelPositions[k][2];
  }
  positions *= eigenVectors.GetVnlMatrix();

  // More CSV data now we know numBoundaryVoxels
  m_FinalStatisticsLog << numBoundaryVoxels << ", " << voxelCount << ", " << minimum[0] << ", ";
  m_FinalStatisticsLog << minimum[1] << ", " << minimum[2] << ", " << maximum[0] << ", ";
  m_FinalStatisticsLog << maximum[1] << ", " << maximum[2] << ", " << mean[0] << ", ";
  m_FinalStatisticsLog << mean[1] << ", " << mean[2] << ", " << boundingBoxDiam << ", ";
  m_FinalStatisticsLog << voxelCountVolume << ", " << voxelCountDia << ", ";

  /* Find max and min of columns from positions matrix */
  itk::Vector<double,3> pmax, pmin, pmean;
  pmax.Fill(INT_MIN);
  pmin.Fill(INT_MAX);
  pmean.Fill(0);
  for (int k = 0; k < voxelCount; k++) {
    if (pmax[0] < positions[k][0]) pmax[0] = positions[k][0];
    if (pmax[1] < positions[k][1]) pmax[1] = positions[k][1];
    if (pmax[2] < positions[k][2]) pmax[2] = positions[k][2];
    if (pmin[0] > positions[k][0]) pmin[0] = positions[k][0];
    if (pmin[1] > positions[k][1]) pmin[1] = positions[k][1];
    if (pmin[2] > positions[k][2]) pmin[2] = positions[k][2];
    pmean += positions[k];
  }    
  pmean /= voxelCount;

  /* Calculate some PCA numbers */
  itk::Vector<double,3> pcaLength = pmax - pmin + 1;

  itk::Matrix<double,2,3> pcaSegment;
  pcaSegment[0][0] = pmin[0];
  pcaSegment[1][0] = pmax[0];
  pcaSegment[0][1] = pcaSegment[1][1] = pmean[1];
  pcaSegment[0][2] = pcaSegment[1][2] = pmean[2];
  pcaSegment *= eigenVectors.GetTranspose();

  /* Calculate convex hull */
  vtkPointsProjectedHull *hull = vtkPointsProjectedHull::New();
  hull->SetDataTypeToDouble();
  hull->SetNumberOfPoints(voxelCount);
  hull->Allocate(sizeof(double),voxelCount);
  for (int k = 0; k < voxelCount; k++) {
    hull->InsertPoint(k,positions[k][0],positions[k][1],positions[k][2]);
  }
  hull->Update();

  // We need the hull points to calculate the hull area
  int numHullPoints = hull->GetSizeCCWHullZ();

  double *pts = new double[2*numHullPoints];
  hull->GetCCWHullZ(pts, numHullPoints);

  /* Set up data structures to calculate area */
  vtkPoints *points = vtkPoints::New();
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(numHullPoints);
  points->Allocate(sizeof(double),numHullPoints);
  for (int k = 0; k < numHullPoints; k++) {
    points->InsertPoint(k,pts[2*k],pts[2*k+1],0.0);
  }
  vtkIdType *ids = new vtkIdType[numHullPoints];
  for (int k = 0; k < numHullPoints; k++) {
    ids[k] = k;
  }

  double normal[3];
  double area = vtkPolygon::ComputeArea(points, numHullPoints,
                                        ids, normal);

  // Write PCA hull area and diameter values in CSV data
  m_FinalStatisticsLog << area*m_AreaFactor << ", " << 
    m_VoxelSize*2.0*sqrt(area/vtkMath::DoublePi()) << ", ";

  delete[] pts;
  delete[] ids;

  /* Ellipsoid fitting - Values not passed to user so comment out for now*/
  double semi_radii[3];
  semi_radii[0] = sqrt(eigenValues[0])*2.0;
  semi_radii[1] = sqrt(eigenValues[1])*2.0;
  semi_radii[2] = sqrt(eigenValues[2])*2.0;
  vtkParametricEllipsoid *ellipsoid = vtkParametricEllipsoid::New();
  ellipsoid->SetXRadius(semi_radii[0]);
  ellipsoid->SetYRadius(semi_radii[1]);
  ellipsoid->SetZRadius(semi_radii[2]);

  ellipsoidVertices.SetSize(numAzimuthal*(numPolar+1),3);

  /* Calculate surface vertices from azimuthal (u) and polar (v) angles
   * Reminder: Azimuthal is clockwise from +ve z axis, 
   * Polar is anti-clockwise from +ve x axis */
  double uvw[3], point[3], Duvw[9];
  uvw[2] = 0;
  for (int iu = 0; iu < numAzimuthal; iu++) {
    uvw[0] = iu*2.0*vtkMath::DoublePi()/double(numAzimuthal);
    for (int ip = 0; ip < numPolar+1; ip++) { // Need ip=numPolar to get lowest point
      uvw[1] = ip*vtkMath::DoublePi()/double(numPolar);
      ellipsoid->Evaluate(uvw, point, Duvw);
      ellipsoidVertices[iu*(numPolar+1)+ip][0] = point[0];
      ellipsoidVertices[iu*(numPolar+1)+ip][1] = point[1];
      ellipsoidVertices[iu*(numPolar+1)+ip][2] = point[2];
    }
  }
  
  ellipsoidVertices *= eigenVectors.GetTranspose();
  for (int k = 0; k < numAzimuthal*(numPolar+1); k++) {
    ellipsoidVertices[k][0] += mean[0];
    ellipsoidVertices[k][1] += mean[1];
    ellipsoidVertices[k][2] += mean[2];
  }

  ellipsoid->Delete();

  points->Delete();
  hull->Delete();

  /* Surface based quantification */
  // Create a Contour Filter 
  vtkSynchronizedTemplates3D *ig = vtkSynchronizedTemplates3D::New();
  /* Set contour value */
  ig->SetValue(0,0.5);
  ig->SetComputeNormals(0);
  ig->SetComputeGradients(0);

  // Setup the input
  vtkImageImport *ii = vtkImageImport::New();
  ii->SetDataExtent(0, m_Dimensions[0] + 1, 0, m_Dimensions[1] + 1, 0, 
    m_Dimensions[2] + 1);
  ii->SetWholeExtent(0, m_Dimensions[0] + 1, 0, m_Dimensions[1] + 1, 0, 
    m_Dimensions[2] + 1);
  ii->SetDataScalarType(m_InputVolumeScalarType);
  ii->SetNumberOfScalarComponents(m_NumComponents);
  ii->SetDataOrigin(m_ImageOrigin[0],
                    m_ImageOrigin[1],
                    m_ImageOrigin[2]);
  ii->SetDataSpacing(m_ImageSpacing[0],
                     m_ImageSpacing[1],
                     m_ImageSpacing[2]);


  // Set up wrapped image ??? ASK SHENG
  IT *inPtr = (IT *)m_ImageData;
  for (int ik = 0; ik < m_Dimensions[2]; ik++) {
    for (int ij = 0; ij < m_Dimensions[1]; ij++) {
      for (int ii = 0; ii < m_Dimensions[0]; ii++) {
        int index = (ik+1)*(m_Dimensions[1]+2)*(m_Dimensions[0]+2)+(ij+1)*
          (m_Dimensions[0]+2)+ii+1;
        // Only work on current intensity level voxels
        if ((int)*inPtr == m_NextValue) m_WrappedImage[index] = m_NextValue;
        else m_WrappedImage[index] = 0;         
        inPtr++;
      }
    }
  }
  /*for (int ii = 0; ii < m_Dimensions[0]*m_Dimensions[1]*m_Dimensions[2]; ii++) {
    if (m_ImageData[ii] == m_NextValue) m_WrappedImage[ii] = m_NextValue;
    else m_WrappedImage[ii] = 0;
  }*/
  ii->SetImportVoidPointer(m_WrappedImage);
  ig->SetInput((vtkDataObject*)(ii->GetOutput()));

  // get the output
  vtkPolyData *od = ig->GetOutput();
  // run the filter
  ig->Update();

  /* Calculate surface area of iso-surface */
  vtkMassProperties *massProperties = vtkMassProperties::New();
  massProperties->SetInput(od);
  massProperties->Update();
  double isoArea = massProperties->GetSurfaceArea();
  double isoVolume = massProperties->GetVolume();

  // Write out iso-surface values to CSV data
  m_FinalStatisticsLog << /*m_AreaFactor**/isoArea << ", " << /*m_VolumeFactor**/isoVolume << ", ";
  m_FinalStatisticsLog << /*m_VoxelSize**/Vol2Dia(isoVolume) << ", ";
  double sphericity = (pow(vtkMath::DoublePi(),1.0/3.0)*
                       pow(6.0*isoVolume,2.0/3.0))/isoArea;
  double aToVRatio = 0.5*Vol2Dia(isoVolume)*isoArea/isoVolume;
  // Write out last CSV data
  m_FinalStatisticsLog << sphericity << ", " << aToVRatio << endl;


  massProperties->Delete();

  

  ii->Delete();
  ig->Delete();

  m_NextValue++;

  return 1;
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
  finalStatisticsLog << m_FinalStatisticsLog.str();
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
double CCPiQuantification3D<IT>::Vol2Dia(double Volume)
{
  return 2.0 * pow((3.0*Volume)/(4.0*vtkMath::DoublePi()),(1.0/3.0));
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

