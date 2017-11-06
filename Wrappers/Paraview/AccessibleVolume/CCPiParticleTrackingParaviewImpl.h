/**
 * CCPi Paraview plugin of particle tracking 
 * The implementation of algorithm is in the CCPi core folder
 *
 * Author : Mr. Srikanth Nagella
 * Date   : 06.06.2014
 */


#ifndef CCPIPARTICLETRACKINGPARAVIEWIMPL_H
#define CCPIPARTICLETRACKINGPARAVIEWIMPL_H

#include "vtkImageAlgorithm.h"
#include <string>
#include <vector>

class VTK_EXPORT CCPiParticleTrackingParaviewImpl : public vtkImageAlgorithm 
{
public:
  static CCPiParticleTrackingParaviewImpl* New();
  vtkTypeMacro(CCPiParticleTrackingParaviewImpl, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  CCPiParticleTrackingParaviewImpl();
  ~CCPiParticleTrackingParaviewImpl();

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);

  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);

private:
  CCPiParticleTrackingParaviewImpl(const CCPiParticleTrackingParaviewImpl&);  // Not implemented.
  void operator=(const CCPiParticleTrackingParaviewImpl&);  // Not implemented.
};

#endif
