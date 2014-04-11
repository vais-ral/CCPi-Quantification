/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file EllipsoidFitting.hpp
This is the header file for class CCPiEllipsoidFitting<IT>
*/

/**
\brief This is a class that can fit ellipsoids to particles from a 
labelled image and save semi-radii and axes to a text file. 

It can be used with an openMP compiler and the CCPiEllipsoidWorker to process the
labelled image in parallel as follows

\code
  CCPiEllipsoidFitting<IT> ef;
  ef.Initialise((void*)info, (void*)pds);
  ef.CreateVoxelIndexList();
  ef.Prepare();
  ef.PrintSummaryData();  
  ef.WriteCSVData("ellipsoid_data.csv");
  int totalVoxels = ef.GetNumVoxelValues();
  #pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < totalVoxels; i++) {

    // Do the real work
    CCPiEllipsoidWorker<IT> *worker = ef.GetWorker(i);
    if (worker != NULL) {
      if (0 == worker->Run()) {
        #pragma omp critical(writefile)
        {
          worker->WriteCSVData("ellipsoid_data.csv");
        }
      }
      delete worker;
    }
  }
\endcode

The default values for voxel size and the area and volume multiplication factors
are 1.0.

*/

#ifndef __ELLIPSOID_FITTING_HPP
#define __ELLIPSOID_FITTING_HPP

#include <cstring>
#include <iostream>
#include <vector>

#include <itkArray.h>
#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkVector.h>

template <class IT> class CCPiEllipsoidWorker;

template <class IT> class CCPiEllipsoidFitting {

  public:

    /** 
     * Default constructor
     * Sets default values for voxel size and the area and volume multiplication
     * factors to 1.0.
     */
    CCPiEllipsoidFitting() : m_ImageData(NULL),
      m_Dimensions(NULL), m_ImageOrigin(NULL), m_ImageSpacing(NULL),
      m_VoxelSize(1.0), m_AreaFactor(1.0), m_VolumeFactor(1.0), m_MinFeatureSize(100) {}

    /**
     * Destructor. Cleans up allocated memory.
     */
    virtual ~CCPiEllipsoidFitting();

    /**
     * Initialise the class with the plugin information and image data passed
     * to a VolView plugin. Pointers are made void* so this header can be 
     * independent of VolView.
     * 
     * \param info Pointer to a vtkVVPluginInfo object
     * \param pds Pointer to a vtkVVProcessDataStruct object
     */
    void Initialise(void *info, void *pds);

    /**
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
    void Initialise(IT *imageData, int *dimensions, 
      int numComponents, int inputScalarType, double *origin, double *spacing, 
      double minValue, double maxValue);

    /**
     * Set the physical voxel size so output has physical size
     * \param voxelSize Physical voxel size
     */
    void SetVoxelSize(double voxelSize);

    /**
     * Create the list of voxel values and list of x,y,z indices for each value.
     * This is a necessary set up step before any real calculation can be 
     * performed
     */
    void CreateVoxelIndexList(void);

    /**
     * Prepare to iterate through the voxel values. Must be called before
     * GetNextWorker().
     */
    void Prepare();

    /**
     * Print summary data of the volume to cout
     */
    void PrintSummaryData(void);

    /**
     * Get the number of different voxel values
     *
     * \return The number of different voxel values.
     */
    int GetNumVoxelValues() {return m_TotalParticles;}

    /**
     * Get the object to do the fitting for object index i
     *
     * \param i The index of the calculation object required
     * \return Pointer to calculation object which the caller should destroy.
     */
    CCPiEllipsoidWorker<IT>* GetWorker(int i);

    /**
     * Get the next object to do the fitting 
     *
     * \return Pointer to calculation object which the caller should destroy.
     */
    CCPiEllipsoidWorker<IT>* GetNextWorker();

    /**
     * Write statistics calculated so far to file.
     *
     * \param filename Name of file to print to.
     */
    void WriteCSVData(std::string filename);

    /// Make the worker a friend so we don't have to transfer a load of data
    friend class CCPiEllipsoidWorker<IT>;

  private:

    /** The raw image data - so m_ImageData[i] is the value of the i^th voxel.
     * This value should not be changed. Use a separate pointer if you want
     * to iterate with a pointer. */
    IT *m_ImageData;
    
    /** Number of voxels in each direction of image */
    int *m_Dimensions;

    /** Number of components at each voxel in the image */
    int m_NumComponents;

    /** Type of scalar data in the image */
    int m_InputVolumeScalarType;

    /** Coordinates of image origin */
    double *m_ImageOrigin;

    /** Spacing of image in x,y,z directions */
    double *m_ImageSpacing;

    /** Minimum raw value from image. */
    double m_MinValue;

    /** Maximum raw value from image. */
    double m_MaxValue;

    /** Physical voxel size */
    double m_VoxelSize;

    /** Area multiplication factors based on physical voxel size */
    double m_AreaFactor;

    /** Volume multiplication factors based on physical voxel size */
    double m_VolumeFactor;

    /** Array of voxel counts for each voxel value. Voxel value is index of
     * this array */
    itk::Array<int> m_VoxelCounts;

    /** Maximum possible number of different voxel values */
    int m_MaxVoxelVals;

    /** Number of different voxel values */
    int m_TotalParticles;

    /** Total number of non-background voxels */
    int m_SumVoxels;

    /** Label must have more than this number of voxels to be included in calculation */
    int m_MinFeatureSize;

    /** Value of next voxel set to be quantified */
    int m_NextValue;

    /** Array to store arrays that will hold the (i,j,k) coordinates of each 
     * voxel with a given intensity value. The index of this array is the 
     * intensity value for those voxels */
    std::vector< std::vector< itk::FixedArray<int,3> > > m_VoxelArray;

    /** Store formatted results to print or display */
    std::ostringstream m_FinalStatisticsLog;

    /**
    \brief Calculate the maximum, minimum and mean positions from a set of (i,j,k)
    coordinates.

    Use itk::Vector for return objects to make subsequent calculations with them
    more straightforward.

    \param Positions Reference to positions we are working with
    \param Minimum Minimum values for (i,j,k) from positions
    \param Maximum Maximum values for (i,j,k) from positions
    \param Mean Mean values for (i,j,k) from positions
    */
    void GetPositionStats(const std::vector< itk::FixedArray<int,3> >& Positions, 
                          itk::Vector<int,3>& Minimum, itk::Vector<int,3>& Maximum, 
                          itk::Vector<double,3>& Mean);

    /**
    \brief Calculate the covariance matrix from the given array of (i,j,k) 
    coordinates

    \param Positions Reference to positions we are working with
    \param CovarianceMatrix The calculated covariance matrix
    */
    void CalculateCovarianceMatrix(
      const std::vector< itk::FixedArray<int,3> >& Positions, 
      itk::Matrix<double,3,3>& CovarianceMatrix);

    /**
     * Work out which is the next voxel value to work on.
     * Updates m_NextValue and returns that value for good measure!
     *
     * \return The voxel value to work on next
     */
    int GetNextVoxelValue();

    /**
     * Write headers for comma separated variable data table to data stream
     */
    void WriteCSVHeader(void);


};

#endif // __ELLIPSOID_FITTING_HPP
