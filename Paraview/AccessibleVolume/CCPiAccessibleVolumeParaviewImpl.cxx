#include "CCPiAccessibleVolumeParaviewImpl.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkImageAlgorithm.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"

#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiAccessibleVolumeITKImpl.h"
#include "CCPiParaviewUserInterface.h"
#include <math.h>

vtkStandardNewMacro(CCPiAccessibleVolumeParaviewImpl);

//----------------------------------------------------------------------------
CCPiAccessibleVolumeParaviewImpl::CCPiAccessibleVolumeParaviewImpl()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
CCPiAccessibleVolumeParaviewImpl::~CCPiAccessibleVolumeParaviewImpl()
{
}

//----------------------------------------------------------------------------
void CCPiAccessibleVolumeParaviewImpl::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


int CCPiAccessibleVolumeParaviewImpl::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector)
{
  vtkImageData *input = vtkImageData::GetData(inputVector[0]);
  vtkImageData *output = vtkImageData::GetData(outputVector);  
  vtkImageData *maskData = vtkImageData::GetData(inputVector[1]);
  bool isTemporaryMaskData = false;

  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());

  int volumeDims[3];
  float voxelSize[3];
  float origin[3];
  origin[0] = output->GetOrigin()[0];
  origin[1] = output->GetOrigin()[1];
  origin[2] = output->GetOrigin()[2];
  voxelSize[0] = output->GetSpacing()[0];
  voxelSize[1] = output->GetSpacing()[1];
  voxelSize[2] = output->GetSpacing()[2];
  volumeDims[0] = output->GetExtent()[1]+1;
  volumeDims[1] = output->GetExtent()[3]+1;
  volumeDims[2] = output->GetExtent()[5]+1;
  vtkWarningMacro(<<origin[0]<<origin[1]<<origin[2]);
  vtkWarningMacro(<<volumeDims[0]<<volumeDims[1]<<volumeDims[2]);
  vtkWarningMacro(<<voxelSize[0]<<voxelSize[1]<<voxelSize[2]);

  if(maskData==NULL) 
  {
	  maskData = vtkImageData::New();
	  isTemporaryMaskData = true;

	  maskData->CopyStructure(input);
	  maskData->AllocateScalars(input->GetScalarType(), input->GetNumberOfScalarComponents());

	  unsigned char* mdata = (unsigned char*)maskData->GetScalarPointer();
	  for(int i=0; i<volumeDims[0]*volumeDims[1]*volumeDims[2];i++) {
		  mdata[i] = 1;
	  }
  }
  CCPiAccessibleVolumeInputImages *imagesInput = new CCPiAccessibleVolumeInputImages(volumeDims,voxelSize,origin, (unsigned char *)input->GetScalarPointer(), (unsigned char *)maskData->GetScalarPointer());
  float sphereDiameterRangeMin=log(MinSphereDiameter), sphereDiameterRangeMax=log(MaxSphereDiameter), imageResolution =ImageResolution;
  CCPiParaviewUserInterface userInterface(this);
  CCPiAccessibleVolumeITKImpl* algoImpl = new CCPiAccessibleVolumeITKImpl(imagesInput, &userInterface, (unsigned char*)output->GetScalarPointer(), sphereDiameterRangeMin, sphereDiameterRangeMax, NumberOfSpheres, imageResolution); 
  algoImpl->Compute();
  std::map<double,double> avResultMap = algoImpl->GetAccessibleVolume();
  for(std::map<double,double>::iterator itr=avResultMap.begin(); itr!=avResultMap.end(); ++itr)
  {
	  	std::ostringstream message;
		message << itr->first <<"    "<< itr->second<< " Spheres";
		userInterface.LogMessage(message.str());
  }
  output->Modified();
  output->GetPointData()->GetScalars()->Modified();
  if(isTemporaryMaskData)
  {
	  maskData->Delete();
  }
  delete imagesInput;
  delete algoImpl;
  return 1;
}


int CCPiAccessibleVolumeParaviewImpl::FillInputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    // the mask input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  return 1;
}