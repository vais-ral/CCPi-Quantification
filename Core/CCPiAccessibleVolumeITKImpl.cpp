#include "CCPiAccessibleVolumeITKImpl.h"

#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkImageRegionConstIterator.h"

#ifdef _WINDOWS 
#define TYPENAME
#else
#define TYPENAME typename
#endif

CCPiAccessibleVolumeITKImpl::CCPiAccessibleVolumeITKImpl(CCPiAccessibleVolumeInputImages *inputImages, CCPiUserApplicationInterface *userInterface,CCPiImageDataUnsignedChar* outputImage, float logMin, float logMax, int numberOfSpheres, float imageResolution)
{	
	InputData = inputImages;
	UserAppInterface = userInterface;
	OutputImage		 = outputImage;
	SphereDiameterRangeMin = logMin;
	SphereDiameterRangeMax = logMax;
	NumberOfSpheres        = numberOfSpheres;
	ImageResolution		   = imageResolution;
	
	isOutputMemoryOwner = false;
}

CCPiAccessibleVolumeITKImpl::~CCPiAccessibleVolumeITKImpl()
{
	if(isOutputMemoryOwner && OutputImage!= NULL)
		delete OutputImage;
}

std::map<double,double> CCPiAccessibleVolumeITKImpl::GetAccessibleVolume()
{
	return AccessibleVolume;
}


DistanceMapImageType::Pointer CCPiAccessibleVolumeITKImpl::GetDistanceMapOfImageWithDanielsson(ImageType::Pointer inputImage)
{
	typedef itk::DanielssonDistanceMapImageFilter< ImageType, DistanceMapImageType >   DistanceMapFilterType;
	DistanceMapFilterType::Pointer distanceMapFilter = DistanceMapFilterType::New();
	//distanceMapFilter->InputIsBinaryOn();
	distanceMapFilter->SetSquaredDistance(false);
	distanceMapFilter->SetUseImageSpacing(false);
//	distanceMapFilter->SetInsideIsPositive(false);

	// Output of import filter is input to distance map filter
	distanceMapFilter->SetInput(inputImage);
	distanceMapFilter->Update();
	return distanceMapFilter->GetOutput();
}

DistanceMapImageType::Pointer CCPiAccessibleVolumeITKImpl::GetDistanceMapOfImageWithMaurer(ImageType::Pointer inputImage)
{
	typedef itk::SignedMaurerDistanceMapImageFilter< ImageType, DistanceMapImageType >   DistanceMapFilterType;
	DistanceMapFilterType::Pointer distanceMapFilter = DistanceMapFilterType::New();
	//distanceMapFilter->InputIsBinaryOn();
	distanceMapFilter->SetSquaredDistance(false);
	distanceMapFilter->SetUseImageSpacing(false);
//	distanceMapFilter->SetInsideIsPositive(false);

	// Output of import filter is input to distance map filter
	distanceMapFilter->SetInput(inputImage);
	distanceMapFilter->Update();
	return distanceMapFilter->GetOutput();
}

DistanceMapImageType::Pointer CCPiAccessibleVolumeITKImpl::ApplyMaskToInputDistanceMapImage(DistanceMapImageType::Pointer inputImage, ImageType::Pointer maskImage)
{
	typedef itk::MaskImageFilter< DistanceMapImageType, ImageType, DistanceMapImageType > MaskFilterType;

	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	// Image to be masked is output from distance map
	maskFilter->SetInput1(inputImage);
	// Mask image is output from mask import filter
	maskFilter->SetInput2(maskImage);
	maskFilter->Update();
	return maskFilter->GetOutput();
}

ImageType::Pointer CCPiAccessibleVolumeITKImpl::BinaryThresholdImage(DistanceMapImageType::Pointer inputImage, float lowerThresholdValue, float upperThresholdValue)
{
		typedef itk::BinaryThresholdImageFilter< DistanceMapImageType, ImageType > ThresholdFilterType;
        ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
        threshold->SetInput(inputImage);
        threshold->SetInsideValue(itk::NumericTraits<unsigned short>::One);
        threshold->SetOutsideValue(itk::NumericTraits<unsigned short>::Zero);
        threshold->SetLowerThreshold(lowerThresholdValue);
        threshold->SetUpperThreshold(upperThresholdValue);	
		threshold->Update();
		return threshold->GetOutput();
}

ImageType::Pointer CCPiAccessibleVolumeITKImpl::SegmentInputImage(DistanceMapImageType::Pointer inputImage, float thresholdValue)
{
	//Make a copy of the input image
	typedef itk::ImageDuplicator< DistanceMapImageType > DistanceMapImageDuplicatorType;
	DistanceMapImageDuplicatorType::Pointer distMapDuplicator = DistanceMapImageDuplicatorType::New();
	distMapDuplicator->SetInputImage(inputImage);
	distMapDuplicator->Update();

	// Typedefs for calculating the volumes
	typedef itk::ConnectedComponentImageFilter< ImageType, ImageType > ConCompFilterType;
	typedef itk::RelabelComponentImageFilter< ImageType, ImageType > RelabelType;


    	ConCompFilterType::Pointer conCompFilter = ConCompFilterType::New();
	RelabelType::Pointer relabel = RelabelType::New();

	ImageType::Pointer binaryImage = BinaryThresholdImage(distMapDuplicator->GetOutput(), thresholdValue, itk::NumericTraits<DistanceMapImageType::PixelType>::max());

	conCompFilter->SetInput(binaryImage);
	conCompFilter->FullyConnectedOn();
	conCompFilter->Update();

	relabel->SetInput( conCompFilter->GetOutput() );
	relabel->Update();


	return relabel->GetOutput();
}

template<typename T, typename O>
T CCPiAccessibleVolumeITKImpl::CreateDuplicateCopyOfImage(T image)
{
	TYPENAME itk::ImageDuplicator< O >::Pointer duplicateImage = itk::ImageDuplicator< O >::New();
	duplicateImage->SetInputImage(image);
	duplicateImage->Update();
	return duplicateImage->GetOutput();
}

ImageType::Pointer CCPiAccessibleVolumeITKImpl::DilateImage(ImageType::Pointer inputImage)
{
        typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StructuringElementType;
        StructuringElementType structuringElement;
        structuringElement.SetRadius(2.5);
        structuringElement.CreateStructuringElement();
        
        typedef itk::GrayscaleDilateImageFilter <ImageType, ImageType, StructuringElementType>
                GrayscaleDilateImageFilterType;

        GrayscaleDilateImageFilterType::Pointer dilateFilter
                = GrayscaleDilateImageFilterType::New();
        dilateFilter->SetInput(inputImage);
        dilateFilter->SetKernel(structuringElement);
        dilateFilter->Update(); 

		return dilateFilter->GetOutput();
}

std::set<int> CCPiAccessibleVolumeITKImpl::GetNonMaskedSegmentLabelsInImage(ImageType::Pointer inputImage)
{
        std::set<int> selectedGroups;

		ImageType::Pointer maskedData = InputData->GetVolumeMaskData();
        itk::ImageRegionConstIterator< ImageType > maskImageIterator(maskedData,
            maskedData->GetRequestedRegion());
        
		itk::ImageRegionConstIterator< ImageType > inputImageIterator(inputImage,
          inputImage->GetRequestedRegion());
        for (inputImageIterator.GoToBegin(), maskImageIterator.GoToBegin(); !inputImageIterator.IsAtEnd(); ++inputImageIterator, ++maskImageIterator) {
            if ( (inputImageIterator.Get() > 0) && (maskImageIterator.Get() == 0) ) selectedGroups.insert(inputImageIterator.Get());
        }

		return selectedGroups;
}

void CCPiAccessibleVolumeITKImpl::BinarizeImageWithSelectedLabels(ImageType::Pointer inputImage, std::set<int> selectedLabels)
{
        itk::ImageRegionIterator< ImageType > labelImageIter(inputImage,
            inputImage->GetRequestedRegion());
        for (labelImageIter.GoToBegin(); !labelImageIter.IsAtEnd(); ++labelImageIter) {
            if ( selectedLabels.find( labelImageIter.Get() ) != selectedLabels.end() ) {
                labelImageIter.Set(1);
            }
            else {
                labelImageIter.Set(0);
            }
        }
		inputImage->Update();
		return;
}

double CCPiAccessibleVolumeITKImpl::ComputeVolumePathAndLabelOutputImage(ImageType::Pointer inputImage, ImageType::Pointer maskImage, ImageType::PixelType* outputImage,int labelId)
{
        itk::ImageRegionConstIterator< ImageType > inputImageIterator( inputImage,
            inputImage->GetRequestedRegion());
        itk::ImageRegionConstIterator< ImageType > maskImageIterator( maskImage,
            inputImage->GetRequestedRegion());
        double volPath = 0;

        int iOut = 0;
        for (inputImageIterator.GoToBegin(), maskImageIterator.GoToBegin(); !inputImageIterator.IsAtEnd(); ++inputImageIterator, ++maskImageIterator, iOut++) {
            if (inputImageIterator.Get()*maskImageIterator.Get() != 0) {
                volPath++;
                outputImage[iOut] = labelId;
            }	
        }	
		return volPath;
}

double CCPiAccessibleVolumeITKImpl::ComputeVolumePathForGivenRadius(DistanceMapImageType::Pointer inputImage, int sphereIndex, float sphereRadius)
{
	ImageType::Pointer segmentedImage = SegmentInputImage(inputImage,sphereRadius);
	ImageType::Pointer copyOfSegementedImage = CreateDuplicateCopyOfImage<ImageType::Pointer, ImageType>(segmentedImage);
	ImageType::Pointer dilatedImage   = DilateImage(segmentedImage);	
	std::set<int>	   nonMaskedLabels = GetNonMaskedSegmentLabelsInImage(dilatedImage);
	BinarizeImageWithSelectedLabels(copyOfSegementedImage, nonMaskedLabels);
	DistanceMapImageType::Pointer distanceMapOfProcessedImage = GetDistanceMapOfImageWithMaurer(copyOfSegementedImage);
	ImageType::Pointer binaryThresholdedImage = BinaryThresholdImage(distanceMapOfProcessedImage,-1*itk::NumericTraits<DistanceMapImageType::PixelType>::max(),sphereRadius);
	double pathVolume = ComputeVolumePathAndLabelOutputImage(binaryThresholdedImage, InputData->GetVolumeMaskData(), OutputImage->GetImage(), sphereIndex+1);
	return pathVolume;
}

void  CCPiAccessibleVolumeITKImpl::Compute()
{
	AllocateMemoryForOutputImageIfNeeded();
	//Calculate DistanceMap
	DistanceMapImageType::Pointer distanceMapOfInputImage = GetDistanceMapOfImageWithMaurer(InputData->GetVolumeData());
	//Mask the Distance map with the input mask image
	DistanceMapImageType::Pointer maskedDistanceMap = ApplyMaskToInputDistanceMapImage(distanceMapOfInputImage, InputData->GetVolumeMaskData());

    // From experience first distance map takes 3 times as long as calculation for
    // one sphere diameter so
    float progressIncrement = 1.0/float(NumberOfSpheres+3);

	if(UserAppInterface->isCancel()) return;
	// For each sphere calculate the path volume fraction and update the AccessibleVolume
	for (unsigned short diameterIdx = 0; diameterIdx < NumberOfSpheres; diameterIdx++) 
	{
		if(UserAppInterface->isCancel()) return;

		std::ostringstream status;
		status << "Processing "<<diameterIdx+1<<" of "<< NumberOfSpheres << " Spheres";
		UserAppInterface->SetProgressValue((3.0+diameterIdx)*progressIncrement);
		UserAppInterface->SetStatusMessage(status.str());
		
		float selectedDiameter = exp(SphereDiameterRangeMin + diameterIdx*(SphereDiameterRangeMax-SphereDiameterRangeMin)/(NumberOfSpheres-1.0));
		float selectedRadius = 0.5*selectedDiameter/ImageResolution;

		double volumePath = ComputeVolumePathForGivenRadius(maskedDistanceMap, diameterIdx, selectedRadius);
		
		double pathVolumeFraction = volumePath/InputData->getScafoldVolume();
		AccessibleVolume.insert(std::pair<double,double>(selectedDiameter, pathVolumeFraction));		
	}
}

void CCPiAccessibleVolumeITKImpl::CopyImage(ImageType::Pointer inputImage, unsigned char* outputImage)
{
	itk::Size<3> size;
	itk::Index<3> index;
	index.Fill(0);
	size[0] = InputData->getDimensions()[0];
		size[1] = InputData->getDimensions()[1];
			size[2] = InputData->getDimensions()[2];
			itk::ImageRegion<3> region(index,size);
	   itk::ImageRegionConstIterator< ImageType > inputImageIterator( inputImage,
            region);

        int iOut = 0;
        for (inputImageIterator.GoToBegin(); !inputImageIterator.IsAtEnd(); ++inputImageIterator, iOut++) {
			outputImage[iOut] = inputImageIterator.Get();
        }	
}

void CCPiAccessibleVolumeITKImpl::CopyImage(DistanceMapImageType::Pointer inputImage, unsigned char* outputImage)
{
	   itk::ImageRegionConstIterator< DistanceMapImageType > inputImageIterator( inputImage,
            inputImage->GetRequestedRegion());
        int iOut = 0;
        for (inputImageIterator.GoToBegin(); !inputImageIterator.IsAtEnd(); ++inputImageIterator, iOut++) {
			outputImage[iOut] = inputImageIterator.Get();
        }	
}

void CCPiAccessibleVolumeITKImpl::SetOutputImage(CCPiImageDataUnsignedChar* outputImage)
{
	if(outputImage!=NULL)
	{
		//check if memory is owned by this class
		if(isOutputMemoryOwner)
			delete OutputImage;
		OutputImage = outputImage;
		isOutputMemoryOwner = false;
	}
}

void CCPiAccessibleVolumeITKImpl::AllocateMemoryForOutputImageIfNeeded()
{
	long size = InputData->getDimensions()[0]*InputData->getDimensions()[1]*InputData->getDimensions()[2];
	long Dims[3];
	Dims[0]=InputData->getDimensions()[0];Dims[1]=InputData->getDimensions()[1];Dims[2]=InputData->getDimensions()[2];
	if(isOutputMemoryOwner||OutputImage!=NULL) {
	}else{
		OutputImage = new CCPiImageDataUnsignedChar(Dims);
		isOutputMemoryOwner = true;
	}
	for(long idx = 0; idx < size;idx++) OutputImage->GetImage()[idx] = 0;
}


CCPiImageDataUnsignedChar* CCPiAccessibleVolumeITKImpl::GetOutputImage()
{
	return OutputImage;
}

CCPiAccessibleVolumeInputImages* CCPiAccessibleVolumeITKImpl::GetInputImages()
{
	return InputData;
}
