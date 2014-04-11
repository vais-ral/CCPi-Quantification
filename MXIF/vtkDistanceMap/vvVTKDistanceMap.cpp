/*=========================================================================

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/VolViewCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/* This plugin performs an integer Mask operation between the first input and
 * the second input. The second input is expected to be the mask. The pixel type
 * of both images must be integer: {char, unsigned char, short, unsigned short,
 * int, unsigned int, long, unsigned long) (not float or double).
 */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkImageEuclideanDistance.h"
#include "vtkImageImport.h"
#include "vtkImageWriter.h"
#include "vtkPointData.h"

#include "vtkVVPluginAPI.h"

template <class IT>
void vvVTKDistanceMapTemplate(vtkVVPluginInfo *info,
                              vtkVVProcessDataStruct *pds, 
                              IT *)
{
  int *dim = info->InputVolumeDimensions;

  vtkImageImport *ii = vtkImageImport::New();
  ii->SetDataExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
  ii->SetWholeExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
  ii->SetDataScalarType(info->InputVolumeScalarType);
  ii->SetNumberOfScalarComponents(
    info->InputVolumeNumberOfComponents);
  ii->SetDataOrigin(info->InputVolumeOrigin[0],
                    info->InputVolumeOrigin[1],
                    info->InputVolumeOrigin[2]);
  ii->SetDataSpacing(info->InputVolumeSpacing[0],
                     info->InputVolumeSpacing[1],
                     info->InputVolumeSpacing[2]);
  ii->SetImportVoidPointer(pds->inData);

  vtkImageEuclideanDistance *distanceMap = vtkImageEuclideanDistance::New();
  distanceMap->SetDimensionality(3);
  distanceMap->SetInput((vtkDataObject*)(ii->GetOutput()));

  //distanceMap->Update();

  //vtkImageCast *castFilter = vtkImageCast::New();
  //castFilter->SetOutputScalarTypeToUnsignedShort();
  //castFilter->SetInput( distanceMap->GetOutput() );
  //castFilter->Update();


  // get the output, would be nice to have VTK write directly 
  // into the output buffer but... VTK is often broken in that regard
  // so we will try, but check afterwards to see if it worked
  //vtkImageData *od = castFilter->GetOutput();
  vtkImageData *od = distanceMap->GetOutput();
  od->UpdateInformation();
  od->SetExtent(0,0,0,0,0,0);
  od->AllocateScalars();
  int size = dim[0] * dim[1] * pds->NumberOfSlicesToProcess *
    info->InputVolumeNumberOfComponents;
  od->SetExtent(0, dim[0] - 1, 0, dim[1] - 1,
                pds->StartSlice,
                pds->StartSlice + pds->NumberOfSlicesToProcess - 1);
  od->GetPointData()->GetScalars()->SetVoidArray(pds->outData,size,1);

  // run the filter
  od->SetUpdateExtent(od->GetExtent());
  info->UpdateProgress(info,(float)0.7,"Calculating distance map");
  od->Update();

  // did VTK not use our memory?
  if (od->GetScalarPointer() != pds->outData)
    {
    memcpy(pds->outData, od->GetScalarPointer(),
           (od->GetPointData()->GetScalars()->GetMaxId() + 1)*
           od->GetPointData()->GetScalars()->GetDataTypeSize());
    }


  /*std::cout << "Writing distance map\n";
  vtkImageWriter *iw = vtkImageWriter::New();
  iw->SetFileDimensionality(3);
  iw->SetFileName("distanceMap.raw");
  iw->SetInput( castFilter->GetOutput() );
  iw->Write();*/

  std::cout << "Processing complete\n";
  info->UpdateProgress(info,(float)1.0,"Processing Complete");

  ii->Delete();
  distanceMap->Delete();
  //iw->Delete();
  //castFilter->Delete();
}


static int ProcessData(void *inf, vtkVVProcessDataStruct *pds)
{

  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  if( info->InputVolumeNumberOfComponents != 1 )
    {
    info->SetProperty( info, VVP_ERROR, "This filter requires a single-component data set as input" ); 
    return -1;
    }

  switch (info->InputVolumeScalarType) {
    vtkTemplateMacro3(vvVTKDistanceMapTemplate, info, pds, 
                      static_cast<VTK_TT *>(0));
  }
  return 0;

}


static int UpdateGUI(void *inf)
{
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  info->OutputVolumeScalarType = VTK_DOUBLE; //UNSIGNED_SHORT /*info->InputVolumeScalarType*/;

  info->OutputVolumeNumberOfComponents = 1;

  for (int i = 0; i < 3; i++) {
    info->OutputVolumeDimensions[i] = info->InputVolumeDimensions[i];
    info->OutputVolumeSpacing[i] = info->InputVolumeSpacing[i];
    info->OutputVolumeOrigin[i] = info->InputVolumeOrigin[i];
  }

  return 1;
}


extern "C" {
  
void VV_PLUGIN_EXPORT vvVTKDistanceMapInit(vtkVVPluginInfo *info)
{
  vvPluginVersionCheck();

  // setup information that never changes
  info->ProcessData = ProcessData;
  info->UpdateGUI   = UpdateGUI;
  info->SetProperty(info, VVP_NAME, "Distance Map (VTK)");
  info->SetProperty(info, VVP_GROUP, "Utility");
  info->SetProperty(info, VVP_TERSE_DOCUMENTATION,
                                    "Distance Map Transform with VTK Distance Map.");
  info->SetProperty(info, VVP_FULL_DOCUMENTATION,
    "This filters computes a Distance map from a binary image using the VTK Distance Map algorithm");

  info->SetProperty(info, VVP_SUPPORTS_IN_PLACE_PROCESSING, "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_PIECES,   "0");
  info->SetProperty(info, VVP_NUMBER_OF_GUI_ITEMS,          "0");
  info->SetProperty(info, VVP_REQUIRED_Z_OVERLAP,           "0");
  info->SetProperty(info, VVP_PER_VOXEL_MEMORY_REQUIRED,    "2");

  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_SERIES_BY_VOLUMES, "0");
  info->SetProperty(info, VVP_PRODUCES_OUTPUT_SERIES, "0");
  info->SetProperty(info, VVP_PRODUCES_PLOTTING_OUTPUT, "0");
}

}
