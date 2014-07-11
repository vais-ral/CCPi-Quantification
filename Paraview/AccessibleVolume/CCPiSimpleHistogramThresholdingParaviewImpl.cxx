#include "CCPiSimpleHistogramThresholdingParaviewImpl.h"

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

#include "CCPiSimpleHistogramThresholdingITKImpl.h"
#include "CCPiParaviewUserInterface.h"
#include <math.h>

vtkStandardNewMacro(CCPiSimpleHistogramThresholdingParaviewImpl);

//----------------------------------------------------------------------------
CCPiSimpleHistogramThresholdingParaviewImpl::CCPiSimpleHistogramThresholdingParaviewImpl()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
CCPiSimpleHistogramThresholdingParaviewImpl::~CCPiSimpleHistogramThresholdingParaviewImpl()
{
}

//----------------------------------------------------------------------------
void CCPiSimpleHistogramThresholdingParaviewImpl::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


int CCPiSimpleHistogramThresholdingParaviewImpl::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector)
{
  vtkImageData *input = vtkImageData::GetData(inputVector[0]);
  vtkImageData *output = vtkImageData::GetData(outputVector->GetInformationObject(0));  

  output->CopyStructure(input);
//output->AllocateScalars(inputVector[0]->GetInformationObject(0));
  output->AllocateScalars(VTK_UNSIGNED_CHAR, input->GetNumberOfScalarComponents());

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

  double range[2];
  input->GetScalarRange(range);
  std::vector<float> peaks;
  switch(input->GetScalarType())
  {
  case VTK_CHAR:
	  peaks = runThresholding((char*)input->GetScalarPointer(), volumeDims, voxelSize, origin, (float)range[0],(float) range[1], (unsigned char*)output->GetScalarPointer());
	  break;
  case VTK_UNSIGNED_CHAR:
	  peaks = runThresholding((unsigned char*)input->GetScalarPointer(), volumeDims, voxelSize, origin, (float)range[0],(float) range[1], (unsigned char*)output->GetScalarPointer());
	  break;
  case VTK_SHORT:
	  peaks = runThresholding((short*)input->GetScalarPointer(), volumeDims, voxelSize, origin, (float)range[0],(float) range[1], (unsigned char*)output->GetScalarPointer());
	  break;
  case VTK_INT:
	  peaks = runThresholding((int*)input->GetScalarPointer(), volumeDims, voxelSize, origin, (float)range[0],(float) range[1], (unsigned char*)output->GetScalarPointer());
	  break;
  case VTK_FLOAT:
	  peaks = runThresholding((float*)input->GetScalarPointer(), volumeDims, voxelSize, origin, (float)range[0],(float) range[1], (unsigned char*)output->GetScalarPointer());
	  break;
  case VTK_DOUBLE:
	  peaks = runThresholding((double*)input->GetScalarPointer(), volumeDims, voxelSize, origin, (float)range[0],(float) range[1], (unsigned char*)output->GetScalarPointer());
	  break;
  default:
	  break;
  }
  output->Modified();
  output->GetPointData()->GetScalars()->Modified();
  return 1;
}

/**
 * Template function to run the ITK k-means image filter
 * @param data The raw data from the input image
 * @param dims The dimensions of the image, (i,j,k)
 * @param voxelSize Size of a voxel in the image
 * @param origin    Origin position of the image
 * @param min		Minimum intensity value in image
 * @param max		Maximum intensity value in image
 * @param output    Raw data for the output image. Set to NULL if no image
 *                  output required.
 * @return The intensity values of first and last peaks.
 */
template <class IT>
std::vector<float> CCPiSimpleHistogramThresholdingParaviewImpl::runThresholding(
    IT *data, const int *dims, const float *voxelSize, const float *origin,
	const float min, const float max,
    unsigned char *output)
{
	long imgDims[3];
	imgDims[0]=dims[0];imgDims[1]=dims[1];imgDims[2]=dims[2];
	CCPiImageData<IT> inputImage(data, imgDims,false); 
	CCPiSimpleHistogramThresholdingITKImpl<IT> shThresholding(&inputImage, dims, voxelSize, origin, min, max);
	shThresholding.Compute();
	CCPiImageDataUnsignedChar outputImage(output, imgDims, false);
	shThresholding.GetOutputImage(&outputImage);
	return shThresholding.GetPeaks();
}

int CCPiSimpleHistogramThresholdingParaviewImpl::FillInputPortInformation(
  int port, vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

int CCPiSimpleHistogramThresholdingParaviewImpl::FillOutputPortInformation(int port, vtkInformation* info)
{
  // now add our info
  if(port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
 
  return 1;
}
