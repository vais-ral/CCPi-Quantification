#include "CCPiParticleTrackingParaviewImpl.h"

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

#include "CCPiParaviewUserInterface.h"
#include <math.h>

vtkStandardNewMacro(CCPiParticleTrackingParaviewImpl);

//----------------------------------------------------------------------------
CCPiParticleTrackingParaviewImpl::CCPiParticleTrackingParaviewImpl()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
CCPiParticleTrackingParaviewImpl::~CCPiParticleTrackingParaviewImpl()
{
}

//----------------------------------------------------------------------------
void CCPiParticleTrackingParaviewImpl::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


int CCPiParticleTrackingParaviewImpl::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector)
{
  vtkImageData *input = vtkImageData::GetData(inputVector[0]);

  return 1;
}



int CCPiParticleTrackingParaviewImpl::FillInputPortInformation(
  int port, vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

int CCPiParticleTrackingParaviewImpl::FillOutputPortInformation(int port, vtkInformation* info)
{
  // now add our info
  if(port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
 
  return 1;
}