/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file Quan3D.hpp
This is the header file for class CCPiQuantification3D<IT>
*/

/**
\brief This is a class that can calculate several characteristics from a 
labelled image.

Use as follows 

\code
  CCPiQuantification3D<IT> quan3D;
  // Initialise from VolView plugin data
  quan3D.Initialise((void*)info, (void*)pds);
  quan3D.CreateVoxelIndexList();
  quan3D.PrepareForQuantification();
  quan3D.PrintData();  
  while (quan3D.NextValueQuantification() != 0) {
    // Report progress or use data from the quantification object
  }
\endcode

The default values for voxel size and the area and volume multiplication factors
are 1.0.

The following characteristics are calculated:
\li Volume by voxel counts
\li Equivalent sphere diameter by voxel counts
\li Bounding box diagonal
\li Principal Component Analysis
\li Ellipsoid fitting by PCA
\li Equivalent circle diameter by PCA
\li Isosurface by marching cube
\li Surface area
\li Surface volume
\li Equivalent sphere diameter from surface volume
\li Sphercity
\li Normalised surface area to volume ratio (Radius*Sa/Vol)
*/

#ifndef __QUAN3D_HPP
#define __QUAN3D_HPP

#include "CCPiDefines.h"
#include <cstring>
#include <iostream>
#include <vector>

#include <itkArray.h>
#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkVector.h>

#include "CCPiLabelQuantificationResult.h"

template <class IT> class CCPiQuantificationWorker;

template <class IT> class CCPI_EXPORT CCPiQuantification3D {

  public:

    /** 
     * Default constructor
     * Sets default values for voxel size and the area and volume multiplication
     * factors to 1.0.
     */
    CCPiQuantification3D() : m_ImageData(NULL), m_WrappedImage(NULL),
      m_Dimensions(NULL), m_ImageOrigin(NULL), m_ImageSpacing(NULL),
      m_VoxelSize(1.0), m_AreaFactor(1.0), m_VolumeFactor(1.0), m_MinFeatureSize(100) {}

    /**
     * Destructor. Cleans up allocated memory.
     */
    virtual ~CCPiQuantification3D();

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
    void Initialise(IT *imageData, const int *dimensions, 
      int numComponents, int inputScalarType, const float *origin, const float *spacing, 
      float minValue, float maxValue);

    /**
     * Set the physical voxel size so output has physical size
     * \param voxelSize Physical voxel size
     */
    void SetVoxelSize(float voxelSize);
    
    /**
     * Set the minimum feature size. Labels with fewer than this number of voxels
     * will not be quantified.
     * \param minFeatureSize Minimum feature size.
     */
    void SetMinFeatureSize(int minFeatureSize) {
        m_MinFeatureSize = minFeatureSize;
    }

    /**
     * Create the list of voxel values and list of x,y,z indices for each value.
     * This is a necessary set up step before any real calculation can be 
     * performed
     */
    void CreateVoxelIndexList(void);

    /**
     * Prepare to iterate through the voxel values. Must be called before
     * NextValueQuantification().
     */
    void PrepareForQuantification();

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
     * Get the object to calculate quantification data for object index i
     *
     * \param i The index of the calculation object required
     * \return Pointer to calculation object which the caller should destroy.
     */
    CCPiQuantificationWorker<IT>* GetWorker(int i);

    /**
     * Get the next object to calculate quantification data 
     *
     * \return Pointer to calculation object which the caller should destroy.
     */
    CCPiQuantificationWorker<IT>* GetNextWorker();

    /**
     * Write statistics calculated so far to file.
     *
     * \param filename Name of file to print to.
     */
    void WriteCSVData(std::string filename);

	/**
	 *
	 * \return Quantification result
	 */
	CCPiLabelQuantificationResult* GetQuantificationResult() { return &QuantificationResult; }
	void SetQuantificationResultByWorker(int labelIndex, std::map<std::string, double> result); 

    /// Make the worker a friend so we don't have to transfer a load of data
    friend class CCPiQuantificationWorker<IT>;

  private:

    /** The raw image data - so inPtr[i] is the value of the i^th voxel.
     * This value should not be changed. Use a separate pointer if you want
     * to iterate with a pointer. */
    IT *m_ImageData;

    /** Image data used for iso-surface calculation */
    IT *m_WrappedImage;

    /** Number of voxels in each direction of image */
    int *m_Dimensions;

    /** Number of components at each voxel in the image */
    int m_NumComponents;

    /** Type of scalar data in the image */
    int m_InputVolumeScalarType;

    /** Coordinates of image origin */
    float *m_ImageOrigin;

    /** Spacing of image in x,y,z directions */
    float *m_ImageSpacing;

    /** Minimum raw value from image. */
    float m_MinValue;

    /** Maximum raw value from image. */
    float m_MaxValue;

    /** Physical voxel size */
    float m_VoxelSize;

    /** Area multiplication factors based on physical voxel size */
    float m_AreaFactor;

    /** Volume multiplication factors based on physical voxel size */
    float m_VolumeFactor;

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

	/** Label Quantification Result */
	CCPiLabelQuantificationResult QuantificationResult;

    /** 
    \brief Calculate sphere diameter from the sphere volume.

    Use the formula \f$ 2*(3v/(4\pi))^{\frac{1}{3}} \f$ where \f$ v \f$ is the 
    sphere volume.

    \param volume The sphere volume
    \return Diameter of the sphere
    */
    inline float Vol2Dia(float Volume);

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

#endif // __QUAN3D_HPP
