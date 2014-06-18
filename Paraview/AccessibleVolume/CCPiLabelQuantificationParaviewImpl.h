/**
 * CCPi Paraview plugin of Label Quantification
 * The implementation of algorithm is in the CCPi core folder
 *
 * Author : Mr. Srikanth Nagella
 * Date   : 13.06.2014
 */


#ifndef CCPILABELQUANTIFICATIONPARAVIEWIMPL_H
#define CCPILABELQUANTIFICATIONPARAVIEWIMPL_H

#include "vtkImageAlgorithm.h"
#include <string>
#include <vtkTable.h>
class VTK_EXPORT CCPiLabelQuantificationParaviewImpl : public vtkImageAlgorithm 
{
public:
  static CCPiLabelQuantificationParaviewImpl* New();
  vtkTypeMacro(CCPiLabelQuantificationParaviewImpl, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkGetMacro(MinFeatureSize, double);
  vtkSetMacro(MinFeatureSize, double);

  vtkGetMacro(VoxelSize, double);
  vtkSetMacro(VoxelSize, double);

protected:
  CCPiLabelQuantificationParaviewImpl();
  ~CCPiLabelQuantificationParaviewImpl();

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);

  double MinFeatureSize;
  double VoxelSize;

  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);

private:
  CCPiLabelQuantificationParaviewImpl(const CCPiLabelQuantificationParaviewImpl&);  // Not implemented.
  void operator=(const CCPiLabelQuantificationParaviewImpl&);  // Not implemented.

template <class IT>
void runQuantification(IT *data, float origin[3],int volumeDims[3], float voxelSize[3],float min,float max, int vtkDataType, vtkTable *outputTable);
};

#endif
