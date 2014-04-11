/** 
 * @file Header file for module using ITK k-means filter.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#ifndef CCPIKMEANSFILTERITK_H
#define CCPIKMEANSFILTERITK_H

#include <vector>

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntSlider.h>

#include "tutorialAPI.h"

class HxUniformScalarField3; // Forward declaration

/** 
 * A simple module that uses the ITK k-means image filter to determine intensity
 * threshold classes. It has user interface to set number of classes.
 */
class TUTORIAL_API CCPiKMeansFilterITK : public HxCompModule
{
    HX_HEADER(CCPiKMeansFilterITK);
    
  public:

    CCPiKMeansFilterITK();
    ~CCPiKMeansFilterITK();

    /** Port providing a button to click to run the module */
    HxPortDoIt portAction;
    
    /** Port providing int text input field with slider for number of classes */
    HxPortIntSlider portClasses;

    /** Perform the calculation in this module. Called by Avizo */
    virtual void compute();
    
  private:

    /**
     * Template function to run the ITK k-means image filter
     * @param data      Raw data from the input image
     * @param dims      Dimensions of the image, (i,j,k)
     * @param voxelSize Size of a voxel in the image
     * @param origin    Origin position of the image
     * @param output    Raw data for the output image. Set to NULL if no image
     *                  output required.
     * @return The estimated means of the threshold classes.
     */
    template<class IT> std::vector<float> runKMeansFilter(IT *data, const int *dims,
                                                          const float *voxelSize, 
                                                          const float *origin,
                                                          unsigned char *output);
  
    /**
     * Create an output with same size as input. The ITK class used for k-means
     * analysis means type must be unsigned char.
     * Check if there is a result that can be re-used. If so check type and size
     * match current input.
     * @param field The input field
     * @return The output to use for showing data in application
     */
    HxUniformScalarField3* createOutput(HxUniformScalarField3 *field);

};

#endif // CCPIKMEANSFILTERITK_H
