#include "CCPiNexusReaderParaviewImpl.h"

#include "vtkStreamingDemandDrivenPipeline.h"
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
#include <math.h>
#include <vtksys/SystemTools.hxx>
#include "NexusWidget/CCPiNexusWidgetDialog.h"
#include "NexusWidget/CCPiNexusReader.h"

vtkStandardNewMacro(CCPiNexusReaderParaviewImpl);


//----------------------------------------------------------------------------
void CCPiNexusReaderParaviewImpl::SetFileName(char *filename)
{
  if (this->FileName == NULL && filename == NULL)
    {
    return;
    }
  if (this->FileName && filename && (!strcmp(this->FileName,filename)))
    {
    return;
    }
  delete [] this->FileName;
  this->FileName = NULL;

  if (filename)
    {
    this->FileName = vtksys::SystemTools::DuplicateString(filename);
    }
//  CCPiNexusWidgetDialog nexusWidget(this->FileName);
//  nexusWidget.exec();
  VariableNameList.clear();
//  for(int i=0;i<nexusWidget.GetSelectedDataSetCount();i++)
//	  VariableNameList.push_back(nexusWidget.GetSelectedDataSet(i));
  CCPiNexusReader reader(this->FileName);
  std::vector<std::string> vecVariableNames = reader.GetVariableNames();
  for(int i=0; i<vecVariableNames.size(); i++)
	  VariableNameList.push_back(vecVariableNames.at(i));
  if(VariableNameList.size()>0)
	this->VariableName = vtksys::SystemTools::DuplicateString(VariableNameList.at(0).c_str());
  this->Modified();
}

//----------------------------------------------------------------------------
CCPiNexusReaderParaviewImpl::CCPiNexusReaderParaviewImpl()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);

  this->FileName                 = NULL;
  this->VariableName			 = NULL;
}

//----------------------------------------------------------------------------
CCPiNexusReaderParaviewImpl::~CCPiNexusReaderParaviewImpl()
{
}

//----------------------------------------------------------------------------
void CCPiNexusReaderParaviewImpl::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


int CCPiNexusReaderParaviewImpl::RequestData(vtkInformation *request,
											 vtkInformationVector **inputVector,
											 vtkInformationVector *outputVector)
{

	if(this->VariableName==NULL||strcmp(this->VariableName,"")==0) return 0;
	CCPiNexusReader reader(this->FileName); 
	int ndims;
	int *dims;
    this->UpdateProgress(0.0);
	ndims = reader.GetDataNumberOfDimensions(this->VariableName);
	dims = new int[ndims];
	reader.GetDataDimensions(std::string(this->VariableName), dims);
	CCPiNexusReader::DATATYPE type = reader.GetDataType(this->VariableName);
	
//	reader.ReadCompleteData(this->VariableName, &data, &ndims, &dims, &dataType,&axisData);
//	vtkGenericWarningMacro(<<"Number of Dimensions"<<ndims);
//	vtkGenericWarningMacro(<<"Size"<<dims[0]<<"x"<<dims[1]<<"x"<<dims[2]);

	vtkImageData *output = vtkImageData::GetData(outputVector->GetInformationObject(0));  
	output->SetOrigin(0,0,0);
	output->SetSpacing(1.0,1.0,1.0);
	output->SetExtent(0, dims[2]-1, 0, dims[1]-1, 0, dims[0]-1);
	output->SetNumberOfScalarComponents(1, outputVector->GetInformationObject(0));
	output->SetScalarType(GetVTKType(type),outputVector->GetInformationObject(0));
	output->AllocateScalars(GetVTKType(type), 1);
	this->UpdateProgress(0.2);
	reader.ReadCompleteDataNoAllocation(this->VariableName, output->GetScalarPointer());
	output->Modified();
    output->GetPointData()->GetScalars()->Modified();
	output->ComputeBounds();
	this->UpdateProgress(1.0);	
	return 1;
}


int CCPiNexusReaderParaviewImpl::FillInputPortInformation(
  int port, vtkInformation* info)
{
	return 1;
}

int CCPiNexusReaderParaviewImpl::FillOutputPortInformation(int port, vtkInformation* info)
{
  // now add our info
	if(port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
 
  return 1;
}

int CCPiNexusReaderParaviewImpl::RequestInformation(  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),  vtkInformationVector *outputVector)
{
	if(this->VariableName==NULL) return 1;
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
	CCPiNexusReader reader(this->FileName); 
	int ndims;
	int *dims;
	ndims = reader.GetDataNumberOfDimensions(this->VariableName);
	dims = new int[ndims];
	reader.GetDataDimensions(std::string(this->VariableName), dims);
	int extent[6];
	extent[0]=extent[2]=extent[4]=0;
	extent[1]=dims[2];extent[3]=dims[1];extent[5]=dims[0];
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
  return 1;
}

int CCPiNexusReaderParaviewImpl::GetNumberOfVariableNameArrays()
{
	return (int)VariableNameList.size();
}

int CCPiNexusReaderParaviewImpl::GetVTKType(CCPiNexusReader::DATATYPE type)
{
	switch(type)
	{
	case CCPiNexusReader::CHAR:
		return VTK_CHAR;
	case CCPiNexusReader::UCHAR:
		return VTK_UNSIGNED_CHAR;
	case CCPiNexusReader::SHORT:
		return VTK_SHORT;
	case CCPiNexusReader::USHORT:
		return VTK_UNSIGNED_SHORT;
	case CCPiNexusReader::INT:
		return VTK_INT;
	case CCPiNexusReader::UINT:
		return VTK_UNSIGNED_INT;
	case CCPiNexusReader::LONG:
		return VTK_LONG;
	case CCPiNexusReader::ULONG:
		return VTK_UNSIGNED_LONG;
	case CCPiNexusReader::FLOAT:
		return VTK_FLOAT;
	case CCPiNexusReader::DOUBLE:
		return VTK_DOUBLE;
	default:
		return VTK_DOUBLE;
	}
}