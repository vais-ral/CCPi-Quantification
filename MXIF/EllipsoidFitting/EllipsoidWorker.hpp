/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file EllipsoidWorker.hpp
This is the header file for class CCPiEllipsoidWorker<IT>
*/

/**
\brief This is a class that can fit an ellipsoid to one labelled region from a 
labelled image and save the semi-radii and axes to a text file. 
It is a friend class to CCPiEllipsoidFitting<IT>
*/

#ifndef __ELLIPSOID_WORKER_HPP
#define __ELLIPSOID_WORKER_HPP

#include "EllipsoidFitting.hpp"

template <class IT> class CCPiEllipsoidWorker {

  public:

    /**
     * Constructor with the CCPiEllipsoidFitting object that holds the data 
     * \param object CCPiEllipsoidFitting object that holds the data we need.
     * \param id The label for the region we are working on.
     */
    CCPiEllipsoidWorker(CCPiEllipsoidFitting<IT> *object, int id) : 
      m_FittingData(object), m_Id(id) {}

    /**
     * Destructor. Cleans up allocated memory.
     */
    virtual ~CCPiEllipsoidWorker();

    /**
     * Do the calculation
     * \return 0 if calculation was run or 1 if this value had insufficient values
     * to be significant.
     */
    int Run();

    /**
     * Write statistics calculated so far to file.
     *
     * \param filename Name of file to print to.
     */
    void WriteCSVData(std::string filename);
    
    /**
     * Draw ellipsoid for this data
     * 
     * \param outputImage Output image raw data array
     */
    void DrawEllipsoid(unsigned short *outputImage);


  private:

    /** Object that holds the data we need */
    CCPiEllipsoidFitting<IT> *m_FittingData;

    /** Id of this worker. Used to obtain data */
    int m_Id;
    
    /** Semi-radii of ellipsoid */
    double m_SemiRadii[3];
    
    /** Axes of ellipsoid in the rows of this matrix */
    double m_Orientation[9];
    
    /** Centroid of ellipsoid */
    int m_Centre[3];

    /** Store formatted results to print or display */
    std::ostringstream m_FinalStatisticsLog;

};

#endif // __ELLIPSOID_WORKER_HPP
