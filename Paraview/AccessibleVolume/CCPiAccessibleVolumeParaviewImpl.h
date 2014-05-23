/**
 * CCPi Paraview plugin of Accessible Volume 
 * The implementation of algorithm is in the CCPi core folder
 *
 * Author : Mr. Srikanth Nagella
 * Date   : 22.05.2014
 */


#ifndef CCPIACCESSIBLEVOLUMEPARAVIEWIMPL_H
#define CCPIACCESSIBLEVOLUMEPARAVIEWIMPL_H

#include "vtkImageAlgorithm.h"
#include <string>

class VTK_EXPORT CCPiAccessibleVolumeParaviewImpl : public vtkImageAlgorithm 
{
public:
  static CCPiAccessibleVolumeParaviewImpl* New();
  vtkTypeMacro(CCPiAccessibleVolumeParaviewImpl, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkGetMacro(MinSphereDiameter, double);
  vtkSetMacro(MinSphereDiameter, double);

  vtkGetMacro(MaxSphereDiameter, double);
  vtkSetMacro(MaxSphereDiameter, double);

  vtkGetMacro(ImageResolution, double);
  vtkSetMacro(ImageResolution, double);

  vtkGetMacro(NumberOfSpheres, int);
  vtkSetMacro(NumberOfSpheres, int);

protected:
  CCPiAccessibleVolumeParaviewImpl();
  ~CCPiAccessibleVolumeParaviewImpl();

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);
  double MinSphereDiameter;
  double MaxSphereDiameter;
  double ImageResolution;
  int    NumberOfSpheres;

  int FillInputPortInformation(int port, vtkInformation* info);

private:
  CCPiAccessibleVolumeParaviewImpl(const CCPiAccessibleVolumeParaviewImpl&);  // Not implemented.
  void operator=(const CCPiAccessibleVolumeParaviewImpl&);  // Not implemented.
};

#endif
