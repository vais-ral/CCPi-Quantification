/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file QuanWorker.cpp
This is the implementation file for class CCPiQuantificationWorker<IT>
*/

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

/**
 * Destructor. Cleans up allocated memory.
 */
template <class IT>
CCPiQuantificationWorker<IT>::~CCPiQuantificationWorker()
{
  m_QuantificationData = NULL;
  if (m_WrappedImage != NULL) delete[] m_WrappedImage;
}


/**
 * Do the calculation
 * \return 0 if calculation was run or 1 if this value had insufficient values
 * to be significant.
 */
template <class IT> 
int CCPiQuantificationWorker<IT>::Run()
{
  int voxelCount;
  double voxelCountVolume, voxelCountDia;
  itk::Vector<int,3> maximum, minimum, boundingBoxSize;
  itk::Vector<double, 3> mean;
  itk::Matrix<double,3,3> covarianceMatrix;
  itk::SymmetricEigenAnalysis< itk::Matrix<double,3,3>, itk::Vector<double,3>,
                               itk::Matrix<double,3,3> > eigenSystem(3);
  itk::Vector<double,3> eigenValues;
  itk::Matrix<double,3,3> eigenVectors;

  itk::VariableSizeMatrix<double> positions;
  int numAzimuthal = 12; /* Azimuthal angle from 0 to 2PI */
  int numPolar = 6; /* Polar angle from 0 to PI */

  itk::VariableSizeMatrix<double> ellipsoidVertices;

  // Create and intialise image to be used in iso-surface construction.
  m_WrappedImage = new IT[(m_QuantificationData->m_Dimensions[0]+2)*
    (m_QuantificationData->m_Dimensions[1]+2)*
    (m_QuantificationData->m_Dimensions[2]+2)];
/*  memset(m_WrappedImage, (IT)0, (m_QuantificationData->m_Dimensions[0]+2)*
    (m_QuantificationData->m_Dimensions[1]+2)*
    (m_QuantificationData->m_Dimensions[2]+2)*sizeof(IT));*/


  /* Do some simple stuff */
  voxelCount = m_QuantificationData->m_VoxelCounts[m_Id];
  // If this feature doesn't have enough voxels to be significant do nothing
  if ((voxelCount == 0) || (voxelCount < m_QuantificationData->m_MinFeatureSize)) {
    return 1;
  }

  voxelCountVolume = voxelCount * m_QuantificationData->m_VolumeFactor;
  voxelCountDia = m_QuantificationData->m_VoxelSize*m_QuantificationData->Vol2Dia(voxelCount);

  // Start writing CSV data
  m_FinalStatisticsLog << m_Id << ", ";
  QuantificationResult["Label"] = m_Id;

  std::vector< itk::FixedArray<int,3> > voxelPositions = m_QuantificationData->m_VoxelArray[m_Id];
  m_QuantificationData->GetPositionStats(voxelPositions, minimum, maximum, mean);

  boundingBoxSize = maximum - minimum + 1;
  double boundingBoxDiam = boundingBoxSize.GetNorm()*m_QuantificationData->m_VoxelSize;
  
  /* Calculate the covariance matrix */  
  m_QuantificationData->CalculateCovarianceMatrix(voxelPositions, covarianceMatrix);

  /* Calculate the eigenvalues and vectors */
  eigenSystem.ComputeEigenValuesAndVectors(covarianceMatrix, eigenValues,
                                           eigenVectors);

  /* Calculate some new positions */
  positions.SetSize(voxelCount,3);
  int numBoundaryVoxels = 0;
  for (int k = 0; k < voxelCount; k++) {
      positions[k][0] = voxelPositions[k][0];
      if ( (voxelPositions[k][0] == 0) || (voxelPositions[k][0] == m_QuantificationData->m_Dimensions[0]-1) ||
           (voxelPositions[k][1] == 0) || (voxelPositions[k][1] == m_QuantificationData->m_Dimensions[1]-1) ||
           (voxelPositions[k][2] == 0) || (voxelPositions[k][2] == m_QuantificationData->m_Dimensions[2]-1) )
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

  QuantificationResult["Border"] = numBoundaryVoxels;
  QuantificationResult["Number_of_Voxels"] = voxelCount;
  QuantificationResult["minX"] = minimum[0];
  QuantificationResult["minY"] = minimum[1];
  QuantificationResult["minZ"] = minimum[2];
  QuantificationResult["maxX"] = maximum[0];
  QuantificationResult["maxY"] = maximum[1];
  QuantificationResult["maxZ"] = maximum[2];
  QuantificationResult["meanX"] = mean[0];
  QuantificationResult["meanY"] = mean[1];
  QuantificationResult["meanZ"] = mean[2];
  QuantificationResult["Bbox_diag"] = boundingBoxDiam;
  QuantificationResult["VoxelCountsVolume"] = voxelCountVolume;
  QuantificationResult["VoxelCountDia"] = voxelCountDia;


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
#if VTK_MAJOR_VERSION > 5 
  m_FinalStatisticsLog << area*m_QuantificationData->m_AreaFactor << ", " << 
    m_QuantificationData->m_VoxelSize*2.0*sqrt(area/vtkMath::Pi()) << ", ";
  QuantificationResult["PCA_2D_area"] = area*m_QuantificationData->m_AreaFactor;
  QuantificationResult["PCA_2D_dia"] =  m_QuantificationData->m_VoxelSize*2.0*sqrt(area/vtkMath::Pi());
#else
  m_FinalStatisticsLog << area*m_QuantificationData->m_AreaFactor << ", " << 
    m_QuantificationData->m_VoxelSize*2.0*sqrt(area/vtkMath::DoublePi()) << ", ";
  QuantificationResult["PCA_2D_area"] = area*m_QuantificationData->m_AreaFactor;
  QuantificationResult["PCA_2D_dia"] =  m_QuantificationData->m_VoxelSize*2.0*sqrt(area/vtkMath::DoublePi());
#endif

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
#if VTK_MAJOR_VERSION > 5
    uvw[0] = iu*2.0*vtkMath::Pi()/double(numAzimuthal);
#else
    uvw[0] = iu*2.0*vtkMath::DoublePi()/double(numAzimuthal);
#endif
    for (int ip = 0; ip < numPolar+1; ip++) { // Need ip=numPolar to get lowest point
#if VTK_MAJOR_VERSION > 5
      uvw[1] = ip*vtkMath::Pi()/double(numPolar);
#else
      uvw[1] = ip*vtkMath::DoublePi()/double(numPolar);
#endif
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

  ii->SetDataExtent(0, m_QuantificationData->m_Dimensions[0] + 1, 
    0, m_QuantificationData->m_Dimensions[1] + 1, 
    0, m_QuantificationData->m_Dimensions[2] + 1);
  ii->SetWholeExtent(0, m_QuantificationData->m_Dimensions[0] + 1, 
    0, m_QuantificationData->m_Dimensions[1] + 1, 
    0, m_QuantificationData->m_Dimensions[2] + 1);
  ii->SetDataScalarType(m_QuantificationData->m_InputVolumeScalarType);
  ii->SetNumberOfScalarComponents(m_QuantificationData->m_NumComponents);
  ii->SetDataOrigin(m_QuantificationData->m_ImageOrigin[0],
                    m_QuantificationData->m_ImageOrigin[1],
                    m_QuantificationData->m_ImageOrigin[2]);
  ii->SetDataSpacing(m_QuantificationData->m_ImageSpacing[0],
                     m_QuantificationData->m_ImageSpacing[1],
                     m_QuantificationData->m_ImageSpacing[2]);


  // Set up wrapped image ??? ASK SHENG
  IT *inPtr = (IT *)m_QuantificationData->m_ImageData;
  for (int ik = 0; ik < m_QuantificationData->m_Dimensions[2]; ik++) {
    for (int ij = 0; ij < m_QuantificationData->m_Dimensions[1]; ij++) {
      for (int ii = 0; ii < m_QuantificationData->m_Dimensions[0]; ii++) {
        int index = (ik+1)*(m_QuantificationData->m_Dimensions[1]+2)*
          (m_QuantificationData->m_Dimensions[0]+2)+(ij+1)*
          (m_QuantificationData->m_Dimensions[0]+2)+ii+1;
        // Only work on current intensity level voxels
        if ((int)*inPtr == m_Id) m_WrappedImage[index] = m_Id;
        else m_WrappedImage[index] = 0;         
        inPtr++;
      }
    }
  }
  ii->SetImportVoidPointer(m_WrappedImage);
  ii->Update();
#if VTK_MAJOR_VERSION > 5 
  ig->SetInputData((vtkDataObject*)ii->GetOutput());
#else
  ig->SetInputConnection((ii->GetOutputPort()));
#endif

  // get the output
  vtkPolyData *od = ig->GetOutput();
  // run the filter
  ig->Update();

  /* Calculate surface area of iso-surface */
  vtkMassProperties *massProperties = vtkMassProperties::New();

#if VTK_MAJOR_VERSION > 5 
  massProperties->SetInputData(ig->GetOutput());
#else
  massProperties->SetInputConnection(ig->GetOutputPort());
#endif

  massProperties->Update();
  double isoArea = massProperties->GetSurfaceArea();
  double isoVolume = massProperties->GetVolume();

  // Write out iso-surface values to CSV data
  m_FinalStatisticsLog << /*m_QuantificationData->m_AreaFactor**/isoArea << 
    ", " << /*m_QuantificationData->m_VolumeFactor**/isoVolume << ", ";
  m_FinalStatisticsLog << /*m_QuantificationData->m_VoxelSize**/m_QuantificationData->Vol2Dia(isoVolume) << ", ";

  QuantificationResult["Surface_Area"] = isoArea;
  QuantificationResult["Surface_Volume"] = isoVolume;
  QuantificationResult["Surface_Dia"] = m_QuantificationData->Vol2Dia(isoVolume);

#if VTK_MAJOR_VERSION > 5 
  double sphericity = (pow(vtkMath::Pi(),1.0/3.0)*
                       pow(6.0*isoVolume,2.0/3.0))/isoArea;
#else
  double sphericity = (pow(vtkMath::DoublePi(),1.0/3.0)*
                       pow(6.0*isoVolume,2.0/3.0))/isoArea;
#endif
  double aToVRatio = 0.5*m_QuantificationData->Vol2Dia(isoVolume)*isoArea/isoVolume;
  // Write out last CSV data
  m_FinalStatisticsLog << sphericity << ", " << aToVRatio << endl;
  QuantificationResult["Sphericity"] = sphericity;
  QuantificationResult["RSaVol"] = aToVRatio;

  massProperties->Delete();

  ii->Delete();
  ig->Delete();


  return 0;

}

/**
 * Write statistics calculated so far to file.
 *
 * \param filename Name of file to print to.
 */
template <class IT>
void CCPiQuantificationWorker<IT>::WriteCSVData(std::string filename)
{
  // Create a file and write to it
  std::ofstream finalStatisticsLog(filename.c_str(), std::ios::app);
  finalStatisticsLog << m_FinalStatisticsLog.str();
  finalStatisticsLog.close();
}


template class CCPiQuantificationWorker<double>;
template class CCPiQuantificationWorker<float>;
template class CCPiQuantificationWorker<long>;
template class CCPiQuantificationWorker<long long>;
template class CCPiQuantificationWorker<unsigned long>;
template class CCPiQuantificationWorker<unsigned long long>;
template class CCPiQuantificationWorker<short>;
template class CCPiQuantificationWorker<unsigned short>;
template class CCPiQuantificationWorker<char>;
template class CCPiQuantificationWorker<unsigned char>;
template class CCPiQuantificationWorker<signed char>;
template class CCPiQuantificationWorker<int>;
template class CCPiQuantificationWorker<unsigned int>;

