#include "CCPiLabelQuantificationParaviewImpl.h"

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

#include "CCPiParaviewUserInterface.h"
#include "CCPiLabelQuantificationResult.h"
#include <math.h>

#include "CCPiLabelQuantificationITKImpl.h"

vtkStandardNewMacro(CCPiLabelQuantificationParaviewImpl);

//----------------------------------------------------------------------------
CCPiLabelQuantificationParaviewImpl::CCPiLabelQuantificationParaviewImpl()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
CCPiLabelQuantificationParaviewImpl::~CCPiLabelQuantificationParaviewImpl()
{
}

//----------------------------------------------------------------------------
void CCPiLabelQuantificationParaviewImpl::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


int CCPiLabelQuantificationParaviewImpl::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector)
{
  vtkImageData *input = vtkImageData::GetData(inputVector[0]);
  vtkTable *outputTable = vtkTable::GetData(outputVector->GetInformationObject(0));

  int volumeDims[3];
  float voxelSize[3];
  float origin[3];
  origin[0] = input->GetOrigin()[0];
  origin[1] = input->GetOrigin()[1];
  origin[2] = input->GetOrigin()[2];
  voxelSize[0] = input->GetSpacing()[0];
  voxelSize[1] = input->GetSpacing()[1];
  voxelSize[2] = input->GetSpacing()[2];
  volumeDims[0] = input->GetExtent()[1]+1;
  volumeDims[1] = input->GetExtent()[3]+1;
  volumeDims[2] = input->GetExtent()[5]+1;

  float min,max;
  min = input->GetScalarRange()[0];
  max = input->GetScalarRange()[1];
  CCPiParaviewUserInterface userInterface(this);

  switch(input->GetScalarType())
  {
	  case VTK_UNSIGNED_CHAR:
		  runQuantification((unsigned char*)input->GetScalarPointer() , origin, volumeDims, voxelSize,  min, max, input->GetScalarType(),outputTable);
		  break;
	  case VTK_UNSIGNED_SHORT:
		  runQuantification((unsigned short*)input->GetScalarPointer() , origin, volumeDims, voxelSize,min, max, input->GetScalarType(),outputTable);
		  break;
	  case VTK_UNSIGNED_INT:
		  runQuantification((unsigned int*)input->GetScalarPointer() , origin, volumeDims, voxelSize,min, max, input->GetScalarType(),outputTable);
		  break;
	  case VTK_INT:
		  runQuantification((int*)input->GetScalarPointer() , origin, volumeDims, voxelSize, min, max, input->GetScalarType(),outputTable);
		  break;
	  default:
		  break;
  }
	  
  return 1;
}


int CCPiLabelQuantificationParaviewImpl::FillInputPortInformation(
  int port, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
	return 1;
}

int CCPiLabelQuantificationParaviewImpl::FillOutputPortInformation(int port, vtkInformation* info)
{
  // now add our info
	if(port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    }
 
  return 1;
}

/**
 * Template function to run calculation
 * @param data          Raw data from the input image
 * @param vtkDataType   Value to indicate IT type to VTK classes.
 *                      Can be found in vtkType.h
 * 
 * @return @TODO
 */
template <class IT>
void CCPiLabelQuantificationParaviewImpl::runQuantification(IT *data, float origin[3], int volumeDims[3], float voxelSize[3], float min,float max, int vtkDataType,vtkTable *output)
{
	long imgDims[3];
	imgDims[0]=volumeDims[0];imgDims[1]=volumeDims[1];imgDims[2]=volumeDims[2];
	CCPiImageData<IT> inputImage(data,imgDims,false);
	CCPiParaviewUserInterface ui(this);
	CCPiLabelQuantificationITKImpl<IT> quantification(&inputImage, &ui,origin,imgDims,voxelSize,min,max,MinFeatureSize,vtkDataType);
	quantification.Compute();
	//Copy the result to output
	CCPiLabelQuantificationResult* result = quantification.GetOutput();
	std::vector<std::string> columnNames = result->GetQuantityNames();
	for(std::vector<std::string>::iterator itr = columnNames.begin(); itr!=columnNames.end();itr++)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName((*itr).c_str());
		output->AddColumn(column);
	}
	std::list<int> labelIndexes = result->GetLabelIndexes();
	output->SetNumberOfRows(labelIndexes.size());
	int columnId=0;
	for(std::vector<std::string>::iterator column_itr = columnNames.begin();column_itr!=columnNames.end();column_itr++, columnId++)
	{
		int rowId=0;
		for(std::list<int>::iterator row_itr=labelIndexes.begin();row_itr!=labelIndexes.end();row_itr++, rowId++)
		{
			double value = result->GetValue(*column_itr, *row_itr);
			output->SetValue(rowId, columnId, vtkVariant(value));
		}
		output->GetColumn(columnId)->Modified();
	}

  output->Modified();

	this->SetProgressText("Processing Complete");
  	this->UpdateProgress(1.0);  
}
