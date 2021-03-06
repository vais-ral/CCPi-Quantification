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
#include <hxcore/internal/HxWorkArea.h>              // Busy-cursor and progress bar
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
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.getSource();
    
    // Get the origin of the data (lower left position of bounding box)
    float origin[3];
    origin[0] = field->getBoundingBox()[0];
	origin[1] = field->getBoundingBox()[2];
	origin[2] = field->getBoundingBox()[4];

	// Get the min and max intensity values in image
	float min = 0.0, max = 0.0;
	field->getRange(min, max);
    
    // Create an output with same size as input. Data type must be unsigned char
    // for ITK filter we are using.
    HxUniformScalarField3 *output = createOutput(field);
    
    // Output shall have same bounding box as input
	output->coords()->setBoundingBox(field->getBoundingBox());
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiSimpleHistogramThresholding", "Computing simple histogram threshold"));

    std::vector<float> peaks;
	int dims[3];
	McDim3l mdims = field->lattice().getDims();
	dims[0] = mdims[0]; dims[1] = mdims[1]; dims[2] = mdims[2];
    switch (field->primType()) {
    
      case McPrimType::MC_UINT8:
        peaks = runThresholding((unsigned char*)field->lattice().dataPtr(),
                                        dims,
                                         field->getVoxelSize().getValue(), origin, 
										 min, max,
                                         (unsigned char*)output->lattice().dataPtr());
        break;
      case McPrimType::MC_INT16:
         peaks = runThresholding((short*)field->lattice().dataPtr(),
                                         dims,
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice().dataPtr());
        break;
      case McPrimType::MC_UINT16:
        peaks = runThresholding((unsigned short*)field->lattice().dataPtr(),
                                         dims,
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice().dataPtr());
        break;
      case McPrimType::MC_INT32:
        peaks = runThresholding((int*)field->lattice().dataPtr(),
                                         dims,
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice().dataPtr());
        break;
      case McPrimType::MC_FLOAT:
        peaks = runThresholding((float*)field->lattice().dataPtr(),
                                         dims,
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice().dataPtr());
        break;
      case McPrimType::MC_DOUBLE:
        peaks = runThresholding((double*)field->lattice().dataPtr(),
                                         dims,
                                         field->getVoxelSize().getValue(), origin, 
										 min, max, 
                                         (unsigned char*)output->lattice().dataPtr());
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
	long imgDims[3];
	imgDims[0]=dims[0];imgDims[1]=dims[1];imgDims[2]=dims[2];
	CCPiImageData<IT> *imgData = new CCPiImageData<IT>(data,imgDims,false);
	CCPiSimpleHistogramThresholdingITKImpl<IT> shThresholding(imgData, dims, voxelSize, origin, min, max);
	shThresholding.Compute();
	CCPiImageDataUnsignedChar *outputImgData = new CCPiImageDataUnsignedChar(output,imgDims,false);
	shThresholding.GetOutputImage(outputImgData); //Copies result image to outputImgData
	delete imgData;
	delete outputImgData;
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
	McDim3l dims = field->lattice().getDims();
    if (output) {
		McDim3l outdims = output->lattice().getDims();
        if ( dims[0] != outdims[0] || dims[1] != outdims[1] ||
            dims[2] != outdims[2] || output->primType() != McPrimType::MC_UINT8 )
            
            output = NULL;
    }
    
    // If necessary create a new result data set
    if (!output) {
        output = new HxUniformScalarField3(dims, McPrimType::MC_UINT8);
        output->composeLabel(field->getName(), "thresholded");
    }
    
    return output;
}