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

#include "omp.h"

#include "Quan3D.hpp"       // Class that control quantification calculation
#include "QuanWorker.hpp"   // Worker class to does calculation for one label

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
  vtkGenericWarningMacro(<<"Min "<<min<<" Max:"<<max);
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
    // Create the controller class for calculations
    CCPiQuantification3D<IT> quan3D;

////    theWorkArea->setProgressValue( 0.01 );
////    theWorkArea->setProgressInfo("Initialising...");


    // Initialise the controller
    quan3D.Initialise(data, volumeDims, 1,
      vtkDataType, origin, voxelSize/****m_field->getVoxelSize().getValue()***/, 
      min, max);

    quan3D.SetMinFeatureSize(/****portMinSize.getValue()*****/MinFeatureSize);
    quan3D.CreateVoxelIndexList();

    quan3D.PrepareForQuantification();

    quan3D.PrintSummaryData();  

////    quan3D.WriteCSVData(portOutputFile.getValue());

    int totalVoxels = quan3D.GetNumVoxelValues(), n = 0;

////    theWorkArea->setProgressValue(0.05);
////    theWorkArea->setProgressInfo("Processing...");

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < totalVoxels; i++) {

        CCPiQuantificationWorker<IT> *worker = NULL;

        // Do the real work
        #pragma omp critical(nextworker)
        {
            worker = quan3D.GetNextWorker();
        }
        if (worker != NULL) {

            if (0 == worker->Run()) {
                #pragma omp critical(writefile)
                {
////                    worker->WriteCSVData(portOutputFile.getValue());
					quan3D.SetQuantificationResultByWorker(worker->GetId(),worker->GetQuantificationResult());
                }
            }
            delete worker;
        }
        #pragma omp atomic
        n++;
        if (omp_get_thread_num() == 0) {
////            theWorkArea->setProgressValue((float)n/(float)totalVoxels);
////            theWorkArea->setProgressInfo("Quantification underway");
        }
    }
////    theMsg->stream() << "Quantification complete" << std::endl;
	//Copy the result to output
	CCPiLabelQuantificationResult* result = quan3D.GetQuantificationResult();
	std::vector<std::string> columnNames = result->GetQuantityNames();
	for(std::vector<std::string>::iterator itr = columnNames.begin(); itr!=columnNames.end();itr++)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName((*itr).c_str());
		output->AddColumn(column);
	}
	std::list<int> labelIndexes = result->GetLabelIndexes();
	output->SetNumberOfRows(labelIndexes.size());
  //outputTable->SetNumberOfRows(avResultMap.size());
  //int rowId=0;
  //for(std::map<double,double>::iterator itr=avResultMap.begin(); itr!=avResultMap.end(); ++itr,rowId++)
  //{
	 // outputTable->SetValue(rowId,0, vtkVariant(itr->first));
	 // outputTable->SetValue(rowId,1, vtkVariant(itr->second));
  //}
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


////    theWorkArea->setProgressValue(1.0);
////    theWorkArea->setProgressInfo("Processing Complete");
    
}
