/*
 *  Template of a compute module
 */

/** 
 * @file Header file for a simple module used for learning about extending Avizo.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#ifndef CCPICOMPUTETHRESHOLD_H
#define CCPICOMPUTETHRESHOLD_H

#include <hxcore/HxCompModule.h>        // Declaration of base class
#include <hxcore/HxPortFloatTextN.h>    // Provides float text input
#include <hxcore/HxPortDoIt.h>          // Provides action button

#include "tutorialAPI.h"    // Export declaration specification

class HxUniformScalarField3; // Forward declaration

/** 
 * A simple threshold module that scans a 3D input image and
 * prints the number of voxels above and below the given threshold values.
 */
class TUTORIAL_API CCPiComputeThreshold : public HxCompModule
{
    // Required macro added by XPand Development Wizard
    HX_HEADER(CCPiComputeThreshold);
    
    /**
     * Create an output with same size and primitive data type as input.
     * Check if there is a result that can be re-used. If so check type, size and
     * datatype match current input.
     * @param field The input field
     * @return The output to use for showing data in application
     */
    HxUniformScalarField3* createOutput(HxUniformScalarField3 *field);

  public:

    CCPiComputeThreshold();
    ~CCPiComputeThreshold();

    /** Port providing a button to click to run the module */
    HxPortDoIt portAction;

    /** Port providing float text input fields */
    HxPortFloatTextN portRange;

    /** Perform the calculation in this module */
    virtual void compute();
};

#endif // CCPICOMPUTETHRESHOLD_H
