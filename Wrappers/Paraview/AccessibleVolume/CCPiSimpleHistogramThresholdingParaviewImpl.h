/**
 * CCPi Paraview plugin of Simple Histogram Thresholding 
 * The implementation of algorithm is in the CCPi core folder
 *
 * Author : Mr. Srikanth Nagella
 * Date   : 06.06.2014
 */


#ifndef CCPISIMPLEHISTOGRAMTHRESHOLDINGPARAVIEWIMPL_H
#define CCPISIMPLEHISTOGRAMTHRESHOLDINGPARAVIEWIMPL_H

#include "vtkImageAlgorithm.h"
#include <string>
#include <vector>

class VTK_EXPORT CCPiSimpleHistogramThresholdingParaviewImpl : public vtkImageAlgorithm 
{
public:
  static CCPiSimpleHistogramThresholdingParaviewImpl* New();
  vtkTypeMacro(CCPiSimpleHistogramThresholdingParaviewImpl, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  CCPiSimpleHistogramThresholdingParaviewImpl();
  ~CCPiSimpleHistogramThresholdingParaviewImpl();

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);

  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);

private:
  CCPiSimpleHistogramThresholdingParaviewImpl(const CCPiSimpleHistogramThresholdingParaviewImpl&);  // Not implemented.
  void operator=(const CCPiSimpleHistogramThresholdingParaviewImpl&);  // Not implemented.
  template <class IT>
  std::vector<float> runThresholding(IT *data, const int *dims, const float *voxelSize, const float *origin,	const float min, const float max,
    unsigned char *output);
};

#endif
