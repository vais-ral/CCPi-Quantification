/** 
 * @file Implementation file for module to calculate accessible volume of masked binary image.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */
#include "api.h"
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
#if   AVIZO_UNSUPPORTED_INTERNAL == 1
#include <hxcore/internal/HxWorkArea.h>
#include <hxrawio/internal/readRawData.h>
#else
#include <hxcore/HxWorkArea.h>
#include <hxrawio/internal/readRawData.h>
#endif

#include <hxfield/HxUniformScalarField3.h>
#if AVIZO_UNSUPPORTED_HXSPREADSHEET == 1
#include <hxspreadsheet/internal/HxSpreadSheet.h>
#else
#include <hxspreadsheet/HxSpreadSheet.h>
#endif

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
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.getSource();
    
    // Check input data
    if (field->primType() != McPrimType::MC_UINT8) {
        theMsg->stream() << "This module only works on a binary image" <<
            " and this input has incorrect data type. Data type must be " <<
            " unsigned char" << std::endl;
        return;
    }
    
    // Get mask data
    HxUniformScalarField3 *maskField = (HxUniformScalarField3*) maskConnection.getSource();
    if (!maskField) {
        theMsg->stream() << "This module requires a second image for masking" <<
            " the input. Please read in masking data and connect to this module" <<
            " in the Project View" << std::endl;
            return;
    }
    
     // Check mask data
    if (maskField->primType() != McPrimType::MC_UINT8) {
        theMsg->stream() << "This module expects a binary mask image" <<
            " and this input has incorrect data type. Data type must be " <<
            " unsigned char" << std::endl;
        return;
    }  
    
    // Get the origin of the data (lower left position of bounding box)
    float origin[3];
    origin[0] = field->getBoundingBox()[0];
	origin[1] = field->getBoundingBox()[2];
	origin[2] = field->getBoundingBox()[4];
    
    // Create an output with same size as input. Data type will be unsigned char
    // as we produce a labelled image.
    HxUniformScalarField3 *output = createOutput(field);
    
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->getBoundingBox());
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorking(
        QApplication::translate("CCPiAccessibleVolume", "Computing accessible volume"));   
        
	//map for volume fraction
	std::map<double,double> outputVolumeFraction;

	int dims[3]; 
	dims[0] = field->lattice().getDims()[0]; dims[1] = field->lattice().getDims()[1]; dims[2] = field->lattice().getDims()[2];
    // Run the main calculation
    run((unsigned char*)field->lattice().dataPtr(),
        dims, field->getVoxelSize().getValue(), origin,
        (unsigned char*)maskField->lattice().dataPtr(), 
        (unsigned char*)output->lattice().dataPtr(),&outputVolumeFraction);
  
    // Register result - adds data object to project view if not already present.
    // Also connects object's master port to compute module
    setResult(output);

	//Register spreadsheet volume fraction
	setResult(createSpreadsheetOutput(field->getName(),outputVolumeFraction));

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
 * @param outputVolumeFraction map of sphere diameter aand volume fraction
 */
void CCPiAccessibleVolume::run(unsigned char *data, const int *dims,
    const float *voxelSize, const float *origin, unsigned char *maskData,
    unsigned char *output,std::map<double,double> *outputVolumeFraction)
{
	long imgDims[3];
	imgDims[0]=dims[0];imgDims[1]=dims[1];imgDims[2]=dims[2];
	CCPiImageDataUnsignedChar* imgData = new CCPiImageDataUnsignedChar(data, imgDims, false);
	CCPiImageDataUnsignedChar* imgMaskData = new CCPiImageDataUnsignedChar(maskData, imgDims, false);
	CCPiAccessibleVolumeInputImages *imgInput = new CCPiAccessibleVolumeInputImages(dims,voxelSize,origin,imgData,imgMaskData);
	CCPiAvizoUserInterface          *userInterface = new CCPiAvizoUserInterface();
	
    // Get data from user interface
    float logMin = log(portDiameterRange.getValue(0));
    float logMax = log(portDiameterRange.getValue(1));
    int numSpheres = portNumSpheres.getValue();

    float imageResolution = portResolution.getValue();

	theWorkArea->setProgressValue( (float)0.05 );

	CCPiImageDataUnsignedChar* imgOutput = new CCPiImageDataUnsignedChar(output, imgDims, false);
	CCPiAccessibleVolumeITKImpl algImpl(imgInput, userInterface, imgOutput,  logMin, logMax, numSpheres, imageResolution);
	algImpl.Compute();

	//createSpreadsheetOutput(algImpl.GetAccessibleVolume());

	writeAccessibleVolumePathFractionToFile(algImpl.GetAccessibleVolume(), portOutputFilename.getFilename().toUtf8().constData());
	//Copy ouput volume fraction 
	std::map<double,double> resultVF = algImpl.GetAccessibleVolume();
	for (std::map<double,double>::iterator resultIterator=resultVF.begin(); resultIterator!=resultVF.end(); ++resultIterator)
	{
		outputVolumeFraction->insert(std::pair<double,double>(resultIterator->first,resultIterator->second));
	}
	delete imgInput;
	delete imgData;
	delete imgMaskData;
	delete imgOutput;
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
        output->composeLabel(field->getName(), "result");
    }
    
    return output;
}

/**
 * Creates the volume path as a spreadsheet and puts in the work area
 * @param name prefix for the module
 * @param volpathMap is map of sphere diameter and its volume path
 * @return the output in a spreadsheet.
 */
HxSpreadSheet* CCPiAccessibleVolume::createSpreadsheetOutput(std::string prefix, std::map<double,double> volpathMap)
{

#if AVIZO_UNSUPPORTED_HXSPREADSHEET == 1
	HxSpreadSheet *output = HxSpreadSheet::createInstance();
#else 
	HxSpreadSheet *output = new HxSpreadSheet();
#endif
	output->addTable(("Accessible Volume("+prefix+")").c_str());
	output->setLabel(("Accessible Volume("+prefix+")").c_str());
	int tableId = 0;//output->findTable(("Accessible Volume("+prefix+")").c_str());
	output->addColumn("Sphere diameter (um)",HxSpreadSheet::Column::FLOAT,tableId);
	output->addColumn("Accessible volume fraction",HxSpreadSheet::Column::FLOAT,tableId);
	int diameterColumnId = output->findColumn("Sphere diameter (um)", HxSpreadSheet::Column::FLOAT,tableId);
	int volumeFractionColumnId = output->findColumn("Accessible volume fraction", HxSpreadSheet::Column::FLOAT,tableId);
	HxSpreadSheet::Column *diameterColumn = output->column(diameterColumnId,tableId);
	HxSpreadSheet::Column *volumeFractionColumn = output->column(volumeFractionColumnId,tableId);
	output->setNumRows(volpathMap.size(),tableId);
	output->setCurrentTableIndex(tableId);
	int rowIndex=0;
	for (std::map<double,double>::iterator resultIterator=volpathMap.begin(); resultIterator!=volpathMap.end(); ++resultIterator, rowIndex++)
	{
		diameterColumn->setValue( rowIndex, resultIterator->first);
		volumeFractionColumn->setValue(rowIndex, resultIterator->second);
	}    
    return output;
}


void CCPiAccessibleVolume::writeAccessibleVolumePathFractionToFile(std::map<double,double> result,std::string fileName)
{
	std::ofstream csvFileWriter;
	csvFileWriter.open(fileName.c_str(),std::ios::trunc);
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
