/**
 * Accessible Volume Input Images is a place holder class for the input images and have utility functions
 * to compute Scafold Volume, Porosity
 *
 * Author: Mr. Srikanth Nagella
 * Date  : 13.05.2014
 */
#ifndef CCPIACCESSIBLEVOLUMEINPUTIMAGES_H
#define CCPIACCESSIBLEVOLUMEINPUTIMAGES_H

#include "itkImage.h"
#include "CCPiImageData.h"

typedef itk::Image< unsigned char, 3 >     ImageType;

class CCPiAccessibleVolumeInputImages
{
public:
	CCPiAccessibleVolumeInputImages(const int *volumeDims, const float *voxelSize, const float *origin, CCPiImageDataUnsignedChar *volumeData, CCPiImageDataUnsignedChar *maskedVolumeData);
	~CCPiAccessibleVolumeInputImages();
	double        getScafoldVolume();
	double		  getScafoldPorosity();
	const unsigned int*	  getDimensions(){ return VolumeDims;}
	const float*  getVoxelSize(){return VoxelSize;}
	const float*  getOrigin(){return Origin;}
	ImageType::Pointer	 GetVolumeMaskData();
	ImageType::Pointer   GetVolumeData();

private:
	CCPiImageDataUnsignedChar *VolumeData;
	CCPiImageDataUnsignedChar *MaskedVolumeData;
	unsigned int VolumeDims[3];
	float		  VoxelSize[3];
	float		  Origin[3];
	double		  ScafoldVolume;
	double		  ScafoldPorosity;

	void calculateScafoldingValues();
	ImageType::Pointer ConvertImageToITKStructure(CCPiImageDataUnsignedChar* data);

};

#endif
