/*
 *  Template of a compute module
 */

#ifndef CCPISIMPLEHISTOGRAMTHRESHOLDING_H
#define CCPISIMPLEHISTOGRAMTHRESHOLDING_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxfield/HxUniformScalarField3.h>

#include "api.h"
#include <vector>

class CCPI_API CCPiSimpleHistogramThresholding : public HxCompModule
{
    HX_HEADER(CCPiSimpleHistogramThresholding);

  public:

//    CCPiSimpleHistogramThresholding();
//    ~CCPiSimpleHistogramThresholding();

    HxPortDoIt portAction;

    virtual void compute();

  private:

    /**
     * Template function to run the thresholding
     * @param data      Raw data from the input image
     * @param dims      Dimensions of the image, (i,j,k)
     * @param voxelSize Size of a voxel in the image
     * @param origin    Origin position of the image
	 * @param min		Minimum intensity value in image
	 * @param max		Maximum intensity value in image
     * @param output    Raw data for the output image. Set to NULL if no image
     *                  output required.
     * @return The intensity values of first and last peaks.
     */
    template<class IT> std::vector<float> runThresholding(IT *data, const int *dims,
                                                          const float *voxelSize, 
                                                          const float *origin,
														  const float min, const float max,
                                                          unsigned char *output);
  
    /**
     * Create an output with same size as input. The type will be unsigned char.
     * Check if there is a result that can be re-used. If so check type and size
     * match current input.
     * @param field The input field
     * @return The output to use for showing data in application
     */
    HxUniformScalarField3* createOutput(HxUniformScalarField3 *field);
};

#endif // CCPISIMPLEHISTOGRAMTHRESHOLDING_H
