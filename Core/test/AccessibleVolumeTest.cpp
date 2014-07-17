/**
 * Testing the Implementation of Accessible Volume using ITK
 * Author : Mr. Srikanth Nagella
 * Date   : 17.07.2014
 */

#include <cmath>
#include "itkImageFileReader.h"
#include "CCPiConsoleUserInterface.h"
#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiAccessibleVolumeITKImpl.h"

int main(int argc, char *argv[])
{
//	Read the input image
	typedef itk::Image<unsigned char, 3> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	
	ReaderType::Pointer DataInputReader = ReaderType::New();
	DataInputReader->SetFileName(argv[1]);
	DataInputReader->GetOutput()->Update();

	ReaderType::Pointer MaskDataInputReader = ReaderType::New();
	MaskDataInputReader->SetFileName(argv[2]);
	MaskDataInputReader->GetOutput()->Update();


	int volumeDims[3];
	volumeDims[0]=volumeDims[1]=volumeDims[2]=128;
	long volumeDimsLong[3];
	volumeDimsLong[0]=volumeDimsLong[1]=volumeDimsLong[2]=128;
	float voxelSize[3];
	voxelSize[0]= voxelSize[1] = voxelSize[2]=1.0;
	float origin[3];
	origin[0]=origin[1]=origin[2]=0;

	CCPiImageDataUnsignedChar DataInput(DataInputReader->GetOutput()->GetBufferPointer(),volumeDimsLong,false);
	CCPiImageDataUnsignedChar MaskDataInput(MaskDataInputReader->GetOutput()->GetBufferPointer(),volumeDimsLong,false);
	CCPiImageDataUnsignedChar Output(volumeDimsLong);
	CCPiAccessibleVolumeInputImages input(volumeDims, voxelSize, origin,&DataInput,&MaskDataInput);

	CCPiConsoleUserInterface console;

	float sphereDiameterRangeMin=log(80.0); 
	float sphereDiameterRangeMax=log(600.0);
	int numberOfSpheres = 11;
	float imageResolution=9.0;

	CCPiAccessibleVolumeITKImpl algo(&input, &console, &Output, sphereDiameterRangeMin, sphereDiameterRangeMax, numberOfSpheres, imageResolution);
	algo.Compute();
	std::map<double,double> outputAV = algo.GetAccessibleVolume();
	if(_isnan(outputAV[80.0])) return 1;
	std::cout<<outputAV[80.0]<<std::endl;
	return 0;
}

