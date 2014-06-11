/**
 * Implementation of Accessible Volume using ITK
 * Author : Mr. Srikanth Nagella
 * Date   : 14.05.2014
 */

#ifndef CCPIACCESSIBLEVOLUMEITKIMPL_H
#define CCPIACCESSIBLEVOLUMEITKIMPL_H

#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiUserApplicationInterface.h"

#include "itkImage.h"


#include <map>
#include <set>

typedef itk::Image< float, 3 >             DistanceMapImageType;

class CCPiAccessibleVolumeITKImpl
{
public:
	CCPiAccessibleVolumeITKImpl(CCPiAccessibleVolumeInputImages *input, CCPiUserApplicationInterface *userAppInterface, unsigned char* outputImage, float sphereDiameterRangeMin, float sphereDiameterRangeMax, int numberOfSpheres, float imageResolution);
	~CCPiAccessibleVolumeITKImpl();
	void                    Compute();
	std::map<double,double> GetAccessibleVolume();
	void	SetOutputImage(unsigned char* outputImage);
	unsigned char*			GetOutputImage();
	CCPiAccessibleVolumeInputImages* GetInputImages();

private:
	CCPiAccessibleVolumeInputImages *InputData;
	CCPiUserApplicationInterface    *UserAppInterface;
	unsigned char*              OutputImage;
	float					    SphereDiameterRangeMin;
	float						SphereDiameterRangeMax;
	int							NumberOfSpheres;
	float						ImageResolution;
	std::map<double,double>     AccessibleVolume;
	bool						isOutputMemoryOwner;

	DistanceMapImageType::Pointer   GetDistanceMapOfImageWithDanielsson(ImageType::Pointer inputImage);
	DistanceMapImageType::Pointer   GetDistanceMapOfImageWithMaurer(ImageType::Pointer inputImage);
	DistanceMapImageType::Pointer	ApplyMaskToInputDistanceMapImage(DistanceMapImageType::Pointer inputImage, ImageType::Pointer maskImage);
	ImageType::Pointer				BinaryThresholdImage(DistanceMapImageType::Pointer inputImage, float lowerThresholdValue, float upperThresholdValue);
	ImageType::Pointer              SegmentInputImage(DistanceMapImageType::Pointer inputImage,float threshold);
	template<typename T,typename O>             
    T        						CreateDuplicateCopyOfImage(T image);
	ImageType::Pointer				DilateImage(ImageType::Pointer inputImage);
	std::set<int>					GetNonMaskedSegmentLabelsInImage(ImageType::Pointer inputImage);
	void							BinarizeImageWithSelectedLabels(ImageType::Pointer inputImage, std::set<int> selectedLabels);
	double							ComputeVolumePathAndLabelOutputImage(ImageType::Pointer inputImage, ImageType::Pointer maskImage, ImageType::PixelType* outputImage,int labelId);
	double							ComputeVolumePathForGivenRadius(DistanceMapImageType::Pointer inputImage, int sphereIndex, float sphereRadius);

	void							CopyImage(ImageType::Pointer inputImage, unsigned char* outputImage);
	void							CopyImage(DistanceMapImageType::Pointer inputImage, unsigned char* outputImage);
	void							AllocateMemoryForOutputImageIfNeeded();
	
};

#endif
