/** 
 * @file Implementation file for module to calculate accessible volume of masked binary image.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#include <fstream>
#include <set>

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxrawio/readRawData.h>

// Classes from ITK
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


#include "CCPiAccessibleVolume.h"

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
    portResolution(this,"resolution",QApplication::translate("CCPiAccessibleVolume", "Image Resolution"),1) 
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
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiAccessibleVolume", "Computing accessible volume"));   
        
    // Run the main calculation
    run((unsigned char*)field->lattice.dataPtr(),
        field->lattice.dims(), field->getVoxelSize().getValue(), origin,
        (unsigned char*)maskField->lattice.dataPtr(), 
        (unsigned char*)output->lattice.dataPtr());
  
    // Register result - adds data object to project view in fot already present.
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
    // Typedefs for image types
    typedef itk::Image< unsigned char, 3 >     ImageType;
    typedef itk::Image< float, 3 >     DistanceMapImageType;
    typedef itk::Image< unsigned char, 3 >     OutputImageType;
    // Typedef for importing image data from VolView data
    typedef itk::ImportImageFilter< unsigned char, 3 > ImportFilterType;

    // Typedef for distance map filters
    typedef itk::DanielssonDistanceMapImageFilter< ImageType, DistanceMapImageType >   FilterType;
    typedef itk::DanielssonDistanceMapImageFilter< OutputImageType, DistanceMapImageType >   SecondFilterType;
    // Typedef for masking the distance filter output
    typedef itk::MaskImageFilter< DistanceMapImageType, ImageType, DistanceMapImageType > MaskFilterType;

    // Typedefs for calculating the volumes
    typedef itk::BinaryThresholdImageFilter< DistanceMapImageType, OutputImageType > ThresholdFilterType;
    typedef itk::ConnectedComponentImageFilter< OutputImageType, OutputImageType > ConCompFilterType;
    typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
    typedef itk::ImageDuplicator< DistanceMapImageType > DistanceMapDuplicatorType;
    typedef itk::ImageDuplicator< OutputImageType > LabelDuplicatorType;

    // Store results as we go along
    std::ostringstream resultsLog;
    
    // Size of image
    int imgSize = dims[0]*dims[1]*dims[2];
    
    int sumROI = 0;
    int sumMaskedImg = 0;

    for (int i = 0; i < imgSize; i++) {
        if ( (int)maskData[i] > 0 ) {
            sumROI++;
            if ( (int)data[i] > 0 ) {
                sumMaskedImg++;
            }
        }
    }
    
    theMsg->stream() << "Volume of scaffold is " << sumROI-sumMaskedImg << std::endl;
    theMsg->stream() << "Scaffold porosity (within ROI) is " << 1.0 - double(sumMaskedImg)/
        double(sumROI) << std::endl;
        
	TYPENAME ImportFilterType::Pointer importFilter = ImportFilterType::New();
    TYPENAME ImportFilterType::Pointer maskImportFilter = ImportFilterType::New();
    TYPENAME FilterType::Pointer distanceMapFilter = FilterType::New();
    //distanceMapFilter->ReleaseDataFlagOn();
    distanceMapFilter->InputIsBinaryOn();

    // Output of import filter is input to distance map filter
    distanceMapFilter->SetInput( importFilter->GetOutput() );
    
    // Data required for importing the image and mask data
    TYPENAME ImportFilterType::IndexType start;
    start.Fill(0);
    TYPENAME ImportFilterType::SizeType size;
    size[0] = dims[0];
    size[1] = dims[1];
    size[2] = dims[2];
    TYPENAME ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(size);
    importFilter->SetSpacing(voxelSize);
    importFilter->SetOrigin(origin);
    importFilter->SetRegion(region);
    maskImportFilter->SetSpacing(voxelSize);
    maskImportFilter->SetOrigin(origin);
    maskImportFilter->SetRegion(region);

    // Set data for full image import filter
    importFilter->SetImportPointer(data, imgSize, false);
    // Set data for mask image import filter
    maskImportFilter->SetImportPointer(maskData, imgSize, false);

    TYPENAME MaskFilterType::Pointer maskFilter = MaskFilterType::New();
    // Image to be masked is output from distance map
    maskFilter->SetInput1(distanceMapFilter->GetOutput());
    // Mask image is output from mask import filter
    maskFilter->SetInput2(maskImportFilter->GetOutput());

    // Output of mask filter is output image
    TYPENAME ImportFilterType::RegionType outputRegion;
    outputRegion.SetIndex(start);
    outputRegion.SetSize(size);
    maskFilter->GetOutput()->SetRegions(outputRegion);

    // Set progress bar, the argument ranges between 0 and 1
    theWorkArea->setProgressValue( (float)0.05 );
    theWorkArea->setProgressInfo("Running initial distance map and mask filter");

    maskFilter->Update();
    
    // Get data from user interface
    float logMin = log(portDiameterRange.getValue(0));
    float logMax = log(portDiameterRange.getValue(1));
    int numSpheres = portNumSpheres.getValue();

    float imageResolution = portResolution.getValue();
    
    // From experience first distance map takes 3 times as long as calculation for
    // one sphere diameter so
    float progressIncrement = 1.0/float(numSpheres+3);

    theWorkArea->setProgressValue(3.0*progressIncrement);

    // Write header for results 
    resultsLog << "Sphere diameter (um), Accessible volume fraction\n";
    
    // Do the volume calulcations
    for (unsigned short diamIdx = 0; diamIdx < numSpheres; diamIdx++) {
        
        float selectedDiameter = exp(logMin + diamIdx*(logMax-logMin)/(numSpheres-1.0));
        float selectedRadius = 0.5*selectedDiameter/imageResolution;
        theMsg->stream() << "Selected diameter and radius " << selectedDiameter <<
            " " << selectedRadius << std::endl;
            
        TYPENAME ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
        TYPENAME ConCompFilterType::Pointer conCompFilter = ConCompFilterType::New();
        TYPENAME RelabelType::Pointer relabel = RelabelType::New();

        theWorkArea->setProgressValue((3.0+diamIdx)*progressIncrement);
        QString message("Working on sphere ");
        QString num;
        num.setNum(diamIdx+1);
        message += num;
        message += " of ";
        num.setNum(numSpheres);
        message += num;
        theWorkArea->setProgressInfo(message);

        // Copy relabelled image out of pipeline for re-use later
        TYPENAME DistanceMapDuplicatorType::Pointer distMapDuplicator = DistanceMapDuplicatorType::New();
        distMapDuplicator->SetInputImage(  maskFilter->GetOutput() );
        distMapDuplicator->Update();

        threshold->SetInput(distMapDuplicator->GetOutput());
        threshold->SetInsideValue(itk::NumericTraits<unsigned short>::One);
        threshold->SetOutsideValue(itk::NumericTraits<unsigned short>::Zero);
        threshold->SetLowerThreshold(selectedRadius);

        conCompFilter->SetInput(threshold->GetOutput());
        conCompFilter->FullyConnectedOn();

        relabel->SetInput( conCompFilter->GetOutput() );
        theWorkArea->setProgressValue((3.0+diamIdx)*progressIncrement);
        relabel->Update();
        theMsg->stream() << "Relabelling complete\n";
        theMsg->stream() << "There are " << relabel->GetNumberOfObjects() << "objects\n";

        // Copy relabelled image out of pipeline for re-use later
        TYPENAME LabelDuplicatorType::Pointer duplicator = LabelDuplicatorType::New();
        duplicator->SetInputImage(  relabel->GetOutput() );
        duplicator->Update();

        typedef itk::BinaryBallStructuringElement<OutputImageType::PixelType, 3> StructuringElementType;
        StructuringElementType structuringElement;
        structuringElement.SetRadius(2.5);
        structuringElement.CreateStructuringElement();
        
        typedef itk::GrayscaleDilateImageFilter <OutputImageType, OutputImageType, StructuringElementType>
                GrayscaleDilateImageFilterType;

        GrayscaleDilateImageFilterType::Pointer dilateFilter
                = GrayscaleDilateImageFilterType::New();
        dilateFilter->SetInput(relabel->GetOutput());
        dilateFilter->SetKernel(structuringElement);
        theWorkArea->setProgressValue((3.0+diamIdx)*progressIncrement);
        dilateFilter->Update(); 

        itk::ImageRegionConstIterator< OutputImageType > iter(dilateFilter->GetOutput(),
          dilateFilter->GetOutput()->GetRequestedRegion());
        
        std::set<int> selectedGroups;
        theWorkArea->setProgressValue((3.0+diamIdx)*progressIncrement);
        itk::ImageRegionConstIterator< ImageType > maskIter(maskImportFilter->GetOutput(),
            maskImportFilter->GetOutput()->GetRequestedRegion());
        for (iter.GoToBegin(), maskIter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++maskIter) {
            if ( (iter.Get() != 0) && (maskIter.Get() == 0) )selectedGroups.insert(iter.Get());
        }
        for (std::set<int>::iterator setIter = selectedGroups.begin(); setIter != selectedGroups.end(); 
             setIter++){
          theMsg->stream() << *setIter << " ";
        }
        theMsg->stream() << std::endl;

        theWorkArea->setProgressValue((3.0+diamIdx)*progressIncrement);
        itk::ImageRegionIterator< OutputImageType > labelImageIter(duplicator->GetOutput(),
            duplicator->GetOutput()->GetRequestedRegion());
        for (labelImageIter.GoToBegin(); !labelImageIter.IsAtEnd(); ++labelImageIter) {
            if ( selectedGroups.find( labelImageIter.Get() ) != selectedGroups.end() ) {
                labelImageIter.Set(1);
            }
            else {
                labelImageIter.Set(0);
            }
        }

        // Values of duplicator output changed above so better call update again  
        duplicator->Update();

        // Do another distance map
        TYPENAME SecondFilterType::Pointer secondDistanceMapFilter = SecondFilterType::New();
        secondDistanceMapFilter->SetInput( duplicator->GetOutput() );
        secondDistanceMapFilter->ReleaseDataFlagOn();
        secondDistanceMapFilter->InputIsBinaryOn(); 

        theWorkArea->setProgressValue((3.0+diamIdx)*progressIncrement);

        TYPENAME ThresholdFilterType::Pointer threshold2 = ThresholdFilterType::New();
        threshold2->SetInput(secondDistanceMapFilter->GetOutput());
        threshold2->SetInsideValue(itk::NumericTraits<unsigned short>::One);
        threshold2->SetOutsideValue(itk::NumericTraits<unsigned short>::Zero);
        threshold2->SetLowerThreshold(0);
        threshold2->SetUpperThreshold(selectedRadius);

        threshold2->Update();
        itk::ImageRegionConstIterator< OutputImageType > iter2(threshold2->GetOutput(),
            threshold2->GetOutput()->GetRequestedRegion());
        int volPath = 0;

        int iOut = 0;
        for (iter2.GoToBegin(), maskIter.GoToBegin(); !iter2.IsAtEnd(); ++iter2, ++maskIter) {
            if (iter2.Get()*maskIter.Get() != 0) {
                volPath++;
                output[iOut] = diamIdx+1;
            }
            iOut++;
        }
        
        theMsg->stream() << "Vol path = " << volPath << std::endl;
        theMsg->stream() << "Path volume fraction = " << double(volPath)/double(sumROI-sumMaskedImg) << std::endl;
        resultsLog << selectedDiameter << ", " << double(volPath)/double(sumROI-sumMaskedImg) << std::endl;
    }

    // Create a results file and write to it
    std::ofstream resultsFile("Accessible_Volume.csv", std::ios::trunc);
    resultsFile << resultsLog.str();
    resultsFile.close();

    theMsg->stream() << "Processing complete\n";
    theWorkArea->setProgressValue((float)1.0);

    std::cout.flush();

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
