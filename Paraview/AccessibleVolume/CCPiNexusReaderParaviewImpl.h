/**
 * This class is a implementation of CCPi Nexus Reader in Paraview
 * Author: Mr. Srikanth Nagella
 * Date  : 29.07.2014
 */

#ifndef CCPINEXUSREADERPARAVIEWIMPL_H
#define CCPINEXUSREADERPARAVIEWIMPL_H

#include "CCPiUserApplicationInterface.h"
#include "vtkImageAlgorithm.h"
#include "NexusWidget/CCPiNexusReader.h"
#include <vector>
#include <string>

class  VTK_EXPORT CCPiNexusReaderParaviewImpl : public vtkImageAlgorithm
{
public:
  static CCPiNexusReaderParaviewImpl* New();
  vtkTypeMacro(CCPiNexusReaderParaviewImpl, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  //Specify file name.
  void SetFileName(char *filename);
  vtkGetStringMacro(FileName);

  //Selected Variable name.
  vtkGetStringMacro(VariableName);
  vtkSetStringMacro(VariableName);

  //Variable Name list
  int GetNumberOfVariableNameArrays();
  const char* GetVariableNameArrayName(int index) { return VariableNameList.at(index).c_str(); }
  int         GetVariableNameArrayStatus(const char* name){return 1;}
  void        SetVariableNameArrayStatus(const char* name, int status){}

protected:
  CCPiNexusReaderParaviewImpl();
  ~CCPiNexusReaderParaviewImpl();

  int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector);

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);

  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);

  char         *FileName;
  char		   *VariableName;
  std::vector<std::string> VariableNameList;

private:
  CCPiNexusReaderParaviewImpl(const CCPiNexusReaderParaviewImpl&);  // Not implemented.
  void operator=(const CCPiNexusReaderParaviewImpl&);  // Not implemented.

  int GetVTKType(CCPiNexusReader::DATATYPE type);
};

#endif