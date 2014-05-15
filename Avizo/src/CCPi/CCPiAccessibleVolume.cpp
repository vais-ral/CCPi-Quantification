/** 
 * @file Implementation file for module to calculate accessible volume of masked binary image.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#include <fstream>
#include <set>

#include <QApplication>

// Classes from ITK
#include "itkIntTypes.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImportImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxrawio/readRawData.h>

#include "CCPiAccessibleVolume.h"

#include "CCPiAccessibleVolumeITKImpl.h"
#include "CCPiAvizoUserInterface.h"

HX_INIT_CLASS(CCPiAccessibleVolume,HxCompModule)

#ifdef _WINDOWS 
  #define TYPENAME
#else
  #define TYPENAME typename
#endif

CCPiAccessibleVolume::CCPiAccessibleVolume() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiAccessibleVolume", "Action")),
    maskConnection(this,"mask",QApplication::translate("CCPiAccessibleVolume", "Mask data")),
    portDiameterRange(this,"diameter",QApplication::translate("CCPiAccessibleVolume", "Min/Max Sphere Diameter (um)"),2),
    portNumSpheres(this,"numSpheres",QApplication::translate("CCPiAccessibleVolume", "Number of Spheres"),1), 
    portResolution(this,"resolution",QApplication::translate("CCPiAccessibleVolume", "Image Resolution"),1),
	portOutputFilename(this,"outputFilename","Output file")
{
    portAction.setLabel(0,"DoIt");
    // Define type of object that can be the masking data
    maskConnection.addType(HxUniformScalarField3::getClassTypeId());
    
    // Define default values for numeric ports
    portDiameterRange.setMinMax(0,1.0,500.0);
    portDiameterRange.setMinMax(1,1.0,1000.0);
    portDiameterRange.setValue(0,80.0);
    portDiameterRange.setValue(1,600.0);
    portNumSpheres.setMinMax(1,20);
    portNumSpheres.setValue(11);
    portResolution.setMinMax(1.0,1000.0);
    portResolution.setValue(9.0);
	portOutputFilename.registerFileType("CSV","csv");
}

CCPiAccessibleVolume::~CCPiAccessibleVolume()
{
}

void CCPiAccessibleVolume::compute()
{
    // Check whether the action port button was clicked
    if (!portAction.wasHit()) return;

    // Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.source();
    
    // Check input data
    if (field->primType() != McPrimType::mc_uint8) {
        theMsg->stream() << "This module only works on a binary image" <<
            " and this input has incorrect data type. Data type must be " <<
            " unsigned char" << std::endl;
        return;
    }
    
    // Get mask data
    HxUniformScalarField3 *maskField = (HxUniformScalarField3*) maskConnection.source();
    if (!maskField) {
        theMsg->stream() << "This module requires a second image for masking" <<
            " the input. Please read in masking data and connect to this module" <<
            " in the Project View" << std::endl;
            return;
    }
    
     // Check mask data
    if (maskField->primType() != McPrimType::mc_uint8) {
        theMsg->stream() << "This module expects a binary mask image" <<
            " and this input has incorrect data type. Data type must be " <<
            " unsigned char" << std::endl;
        return;
    }  
    
    // Get the origin of the data (lower left position of bounding box)
    float origin[3];
    origin[0] = field->bbox()[0];
    origin[1] = field->bbox()[2];
    origin[2] = field->bbox()[4];
    
    // Create an output with same size as input. Data type will be unsigned char
    // as we produce a labelled image.
    HxUniformScalarField3 *output = createOutput(field);
    
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->bbox());
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorking(
        QApplication::translate("CCPiAccessibleVolume", "Computing accessible volume"));   
        
    // Run the main calculation
    run((unsigned char*)field->lattice.dataPtr(),
        field->lattice.dims(), field->getVoxelSize().getValue(), origin,
        (unsigned char*)maskField->lattice.dataPtr(), 
        (unsigned char*)output->lattice.dataPtr());
  
    // Register result - adds data object to project view if not already present.
    // Also connects object's master port to compute module
    setResult(output);

    // Stop progress bar
    theWorkArea->stopWorking();

}

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
void CCPiAccessibleVolume::run(unsigned char *data, const int *dims,
    const float *voxelSize, const float *origin, unsigned char *maskData,
    unsigned char *output)
{
	CCPiAccessibleVolumeInputImages *imgInput = new CCPiAccessibleVolumeInputImages(dims,voxelSize,origin,data,maskData);
	CCPiAvizoUserInterface          *userInterface = new CCPiAvizoUserInterface();

    // Get data from user interface
    float logMin = log(portDiameterRange.getValue(0));
    float logMax = log(portDiameterRange.getValue(1));
    int numSpheres = portNumSpheres.getValue();

    float imageResolution = portResolution.getValue();

	theWorkArea->setProgressValue( (float)0.05 );

	CCPiAccessibleVolumeITKImpl algImpl(imgInput, userInterface, output,  logMin, logMax, numSpheres, imageResolution);
	algImpl.Compute();

	writeAccessibleVolumePathFractionToFile(algImpl.GetAccessibleVolume(), portOutputFilename.getFilename());
	delete imgInput;
	delete userInterface;
}

/**
 * Create an output with same size as input. The output image is essentially
 * a labelled image so data type is unsigned char.
 * Check if there is a result that can be re-used. If so check type and size
 * match current input.
 * @param field The input field
 * @return The output to use for showing data in application
 */
HxUniformScalarField3* CCPiAccessibleVolume::createOutput(HxUniformScalarField3 *field)
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
        output->composeLabel(field->getName(), "result");
    }
    
    return output;
}


void CCPiAccessibleVolume::writeAccessibleVolumePathFractionToFile(std::map<double,double> result,std::string fileName)
{
	std::ofstream csvFileWriter;
	csvFileWriter.open(fileName,std::ios::trunc);
	if(!csvFileWriter.is_open())
	{
		theMsg->stream() << "Error Opening output file"<<std::endl;
	}
	//Header
	csvFileWriter << "Sphere diameter (um), Accessible volume fraction"<<std::endl;
	for (std::map<double,double>::iterator resultIterator=result.begin(); resultIterator!=result.end(); ++resultIterator)
	{
		csvFileWriter << resultIterator->first<<","<<resultIterator->second<<std::endl;
	}
	csvFileWriter.close();
}