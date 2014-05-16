/*
 *  Template of a compute module
 */

#include <QApplication>

// Classes from ITK
#include "itkBinaryThresholdImageFilter.h"
#include "itkImage.h"
#include "itkImportImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>              // Busy-cursor and progress bar
#include <hxfield/HxUniformScalarField3.h>

#include "CCPiSimpleHistogramThresholding.h"
#include "CCPiSimpleHistogramThresholdingITKImpl.h"

HX_INIT_CLASS(CCPiSimpleHistogramThresholding,HxCompModule)

CCPiSimpleHistogramThresholding::CCPiSimpleHistogramThresholding() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiSimpleHistogramThresholding", "Action"))
{
    portAction.setLabel(0,"DoIt");
}

CCPiSimpleHistogramThresholding::~CCPiSimpleHistogramThresholding()
{
}

void CCPiSimpleHistogramThresholding::compute()
{
    // Check whether the action port button was clicked
    if (!portAction.wasHit()) return;

    // Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.source();
    
    // Get the origin of the data (lower left position of bounding box)
    float origin[3];
    origin[0] = field->bbox()[0];
    origin[1] = field->bbox()[2];
    origin[2] = field->bbox()[4];

	// Get the min and max intensity values in image
	float min = 0.0, max = 0.0;
	field->getRange(min, max);
    
    // Create an output with same size as input. Data type must be unsigned char
    // for ITK filter we are using.
    HxUniformScalarField3 *output = createOutput(field);
    
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->bbox());
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiSimpleHistogramThresholding", "Computing simple histogram threshold"));

    std::vector<float> peaks;
    switch (field->primType()) {
    
      case McPrimType::mc_uint8:
        peaks = runThresholding((unsigned char*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
										 min, max,
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_int16:
         peaks = runThresholding((short*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_uint16:
        peaks = runThresholding((unsigned short*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_int32:
        peaks = runThresholding((int*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_float:
        peaks = runThresholding((float*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_double:
        peaks = runThresholding((double*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice.dataPtr());
         break;
      default:
        theMsg->stream() << "The input data has a datatype that this module" <<
                            " cannot handle" << std::endl;
   }
    
    // Register result - adds data object to project view in fot already present.
    // Also connects object's master port to compute module
    setResult(output);

    // Stop progress bar
    theWorkArea->stopWorking();

    // Finally print the result
    theMsg->stream() <<  "Historgram peaks used are " << peaks[0] << " " <<
		peaks[1] << std::endl;

}

/**
 * Template function to run the ITK k-means image filter
 * @param data The raw data from the input image
 * @param dims The dimensions of the image, (i,j,k)
 * @param voxelSize Size of a voxel in the image
 * @param origin    Origin position of the image
 * @param min		Minimum intensity value in image
 * @param max		Maximum intensity value in image
 * @param output    Raw data for the output image. Set to NULL if no image
 *                  output required.
 * @return The intensity values of first and last peaks.
 */
template <class IT>
std::vector<float> CCPiSimpleHistogramThresholding::runThresholding(
    IT *data, const int *dims, const float *voxelSize, const float *origin,
	const float min, const float max,
    unsigned char *output)
{

	CCPiSimpleHistogramThresholdingITKImpl<IT> shThresholding(data, dims, voxelSize, origin, min, max);
	shThresholding.Compute();
	CCPiSimpleHistogramThresholdingITKImpl<IT>::OutputImageType::Pointer image = shThresholding.GetOutputImage();

	return shThresholding.GetPeaks();
}

/**
 * Create an output with same size as input. The type will be unsigned char.
 * Check if there is a result that can be re-used. If so check type and size
 * match current input.
 * @param field The input field
 * @return The output to use for showing data in application
 */
HxUniformScalarField3* CCPiSimpleHistogramThresholding::createOutput(HxUniformScalarField3 *field)
{
    // Check if there is a result which we can reuse
    HxUniformScalarField3 *output = (HxUniformScalarField3*) getResult();
    
    // Check for proper type
    if ( output && !output->isOfType(HxUniformScalarField3::getClassTypeId()) )
        output = NULL;
    
    // Check if size and primitive type still match current input
    const int *dims = field->lattice.dims();
    if (output) {
        const int *outdims = output->lattice.dims();
        if ( dims[0] != outdims[0] || dims[1] != outdims[1] ||
            dims[2] != outdims[2] || output->primType() != McPrimType::mc_uint8 )
            
            output = NULL;
    }
    
    // If necessary create a new result data set
    if (!output) {
        output = new HxUniformScalarField3(dims, McPrimType::mc_uint8);
        output->composeLabel(field->getName(), "thresholded");
    }
    
    return output;
}