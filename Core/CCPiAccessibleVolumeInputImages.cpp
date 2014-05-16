#include "CCPiAccessibleVolumeInputImages.h"

#include "itkImportImageFilter.h"


#ifdef _WINDOWS 
  #define TYPENAME
#else
  #define TYPENAME typename
#endif

CCPiAccessibleVolumeInputImages::CCPiAccessibleVolumeInputImages(const int *volumeDims, const float *voxelSize, const float *origin, unsigned char *volumeData, unsigned char *maskedVolumeData)
{
	VolumeDims[0] = volumeDims[0]; VolumeDims[1] = volumeDims[1]; VolumeDims[2] = volumeDims[2];
	VoxelSize [0] = voxelSize[0];  VoxelSize[1]  = voxelSize[1];  VoxelSize[2]  = voxelSize[2];
	Origin[0]	  = origin[0];     Origin[1]     = origin[1];     Origin[2]     = origin[2];
	VolumeData    = volumeData;
	MaskedVolumeData = maskedVolumeData;
	calculateScafoldingValues();
}

CCPiAccessibleVolumeInputImages::~CCPiAccessibleVolumeInputImages()
{
	VolumeData = 0;
	MaskedVolumeData = 0;
}

void CCPiAccessibleVolumeInputImages::calculateScafoldingValues()
{   
    unsigned long imgSize = VolumeDims[0]*VolumeDims[1]*VolumeDims[2];
    
    unsigned long sumROI = 0;
    unsigned long sumMaskedImg = 0;

    for (unsigned long int i = 0; i < imgSize; i++) {
        if ( (int)MaskedVolumeData[i] > 0 ) {
            sumROI++;
            if ( (int)VolumeData[i] > 0 ) {
                sumMaskedImg++;
            }
        }
    }
	this->ScafoldVolume = sumROI-sumMaskedImg;
	this->ScafoldPorosity = 1 - double(sumMaskedImg)/double(sumROI);
}

double CCPiAccessibleVolumeInputImages::getScafoldVolume()
{
	return this->ScafoldVolume;
}
double CCPiAccessibleVolumeInputImages::getScafoldPorosity()
{
	return this->ScafoldPorosity;
}

ImageType::Pointer CCPiAccessibleVolumeInputImages::GetVolumeMaskData()
{
	return ConvertImageToITKStructure(MaskedVolumeData);
}

ImageType::Pointer CCPiAccessibleVolumeInputImages::GetVolumeData()
{
	return ConvertImageToITKStructure(VolumeData);
}

ImageType::Pointer CCPiAccessibleVolumeInputImages::ConvertImageToITKStructure(unsigned char* data)
{
	typedef itk::ImportImageFilter< unsigned char, 3 > ImportFilterType;

    // Size of image
    unsigned long imgSize = ((unsigned long)VolumeDims[0])*VolumeDims[1]*VolumeDims[2];
	TYPENAME ImportFilterType::Pointer importFilter = ImportFilterType::New();
    TYPENAME ImportFilterType::IndexType start;
    start.Fill(0);
    TYPENAME ImportFilterType::SizeType size;
    size[0] = VolumeDims[0];
    size[1] = VolumeDims[1];
    size[2] = VolumeDims[2];
    TYPENAME ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(size);
    importFilter->SetSpacing(VoxelSize);
    importFilter->SetOrigin(Origin);
    importFilter->SetRegion(region);
	importFilter->SetImportPointer(data, imgSize, false);
	importFilter->Update();
	return importFilter->GetOutput();
}