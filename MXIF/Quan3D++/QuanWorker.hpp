/*=========================================================================

  Copyright (c) David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file QuanWorker.hpp
This is the header file for class CCPiQuantificationWorker<IT>
*/

/**
\brief This is a class that can calculate several characteristics from a 
labelled image. It is a friend class to CCPiQuantification<IT>
*/

#ifndef __QUANWORKER_HPP
#define __QUANWORKER_HPP

#include "Quan3D.hpp"

template <class IT> class CCPiQuantificationWorker {

  public:

    /**
     * Constructor with the CCPiQuantification object that holds the data 
     * \param object CCPiQuantification object that holds the data we need.
     */
    CCPiQuantificationWorker(CCPiQuantification3D<IT> *object, int id) : 
      m_QuantificationData(object), m_Id(id), m_WrappedImage(NULL) {}

    /**
     * Destructor. Cleans up allocated memory.
     */
    virtual ~CCPiQuantificationWorker();

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


  private:

    /** Object that holds the data we need */
    CCPiQuantification3D<IT> *m_QuantificationData;

    /** Id of this worker. Used to obtain data */
    int m_Id;

    /** Image data used for iso-surface calculation */
    IT *m_WrappedImage;

    /** Store formatted results to print or display */
    std::ostringstream m_FinalStatisticsLog;

};

#endif // __QUANWORKER_HPP
