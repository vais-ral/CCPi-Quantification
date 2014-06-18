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
#include "vtkDoubleArray.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"

#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiAccessibleVolumeITKImpl.h"
#include "CCPiParaviewUserInterface.h"
#include <math.h>

vtkStandardNewMacro(CCPiAccessibleVolumeParaviewImpl);

//----------------------------------------------------------------------------
CCPiAccessibleVolumeParaviewImpl::CCPiAccessibleVolumeParaviewImpl()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
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
  vtkImageData *output = vtkImageData::GetData(outputVector->GetInformationObject(0));  
  vtkImageData *maskData = vtkImageData::GetData(inputVector[1]);
  vtkTable *outputTable = vtkTable::GetData(outputVector->GetInformationObject(1));

  bool isTemporaryMaskData = false;

  output->CopyStructure(input);
  output->AllocateScalars(inputVector[0]->GetInformationObject(0));

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
  this->SetProgress(0.5);
  CCPiAccessibleVolumeITKImpl* algoImpl = new CCPiAccessibleVolumeITKImpl(imagesInput, &userInterface, (unsigned char*)output->GetScalarPointer(), sphereDiameterRangeMin, sphereDiameterRangeMax, NumberOfSpheres, imageResolution); 
  algoImpl->Compute();
  std::map<double,double> avResultMap = algoImpl->GetAccessibleVolume();
  vtkSmartPointer<vtkDoubleArray> sphereDiaColumn = vtkSmartPointer<vtkDoubleArray>::New();
  sphereDiaColumn->SetName("Sphere Diameter(um)");
  vtkSmartPointer<vtkDoubleArray> volumeFractionColumn = vtkSmartPointer<vtkDoubleArray>::New();
  volumeFractionColumn->SetName("Accessible Volume Fraction");
  outputTable->AddColumn(sphereDiaColumn);
  outputTable->AddColumn(volumeFractionColumn);
  outputTable->SetNumberOfRows(avResultMap.size());
  int rowId=0;
  for(std::map<double,double>::iterator itr=avResultMap.begin(); itr!=avResultMap.end(); ++itr,rowId++)
  {
	  outputTable->SetValue(rowId,0, vtkVariant(itr->first));
	  outputTable->SetValue(rowId,1, vtkVariant(itr->second));
  }
  outputTable->Modified();
  outputTable->GetColumn(0)->Modified();
  outputTable->GetColumn(1)->Modified();
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

int CCPiAccessibleVolumeParaviewImpl::FillOutputPortInformation(int port, vtkInformation* info)
{
  // now add our info
  if(port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
  else if(port == 1)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    }
 
  return 1;
}