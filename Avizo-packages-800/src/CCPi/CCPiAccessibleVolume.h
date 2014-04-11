/** 
 * @file Header file for module to calculate accessible volume of masked binary image.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#ifndef CCPIACCESSIBLEVOLUME_H
#define CCPIACCESSIBLEVOLUME_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxConnection.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatTextN.h> 
#include <hxcore/HxPortIntTextN.h> 

#include "api.h"

/** 
 * A module that uses various ITK classes to calculate the accessible volume
 * for a binary image with associated mask image.
 * It has user interface to set mask image, range of sphere diameters, number
 * of spheres and the image resolution.
 * 
 * The algorithm is described in detail in ...
 */
class CCPI_API CCPiAccessibleVolume : public HxCompModule
{
    HX_HEADER(CCPiAccessibleVolume);

  public:

    CCPiAccessibleVolume();
    ~CCPiAccessibleVolume();

    /** Port providing a button to click to run the module */
    HxPortDoIt portAction;
    /** Connection to masking data */
    HxConnection maskConnection;
    /** Port providing float text input fields for min and max sphere diameter */
    HxPortFloatTextN portDiameterRange;
    /** Port providing float text input fields for number of sphere diameters */
    HxPortIntTextN portNumSpheres;
    /** Port providing float text input fields for image resolution */
    HxPortFloatTextN portResolution;

    /** Perform the calculation in this module. Called by Avizo. */
    virtual void compute();
    
  private:

    /**
     * Function to run the calculation
     * @param data      Raw data from the input image
     * @param dims      Dimensions of the image, (i,j,k)
     * @param voxelSize Size of a voxel in the image
     * @param origin    Origin position of the image
     * @param maskData  Raw data from the masking image
     * @param output    Raw data for the output image. Set to NULL if no image
     *                  output required.
     */
    void run(unsigned char *data, const int *dims,
             const float *voxelSize, const float *origin,
             unsigned char *maskData, unsigned char *output);
             
    /**
     * Create an output with same size as input. The output image is essentially
     * a labelled image so data type is unsigned char.
     * Check if there is a result that can be re-used. If so check type and size 
     * match current input.
     * @param field The input field
     * @return The output to use for showing data in application
     */
    HxUniformScalarField3* createOutput(HxUniformScalarField3 *field);
};

#endif // CCPIACCESSIBLEVOLUME_H
