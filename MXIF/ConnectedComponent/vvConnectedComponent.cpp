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

#include "itkConnectedComponentImageFilter.h"
#include "itkImage.h"
#include "itkImportImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "vtkVVPluginAPI.h"

template <class IT>
void vvConnectedComponentTemplate(vtkVVPluginInfo *info,
                                vtkVVProcessDataStruct *pds, 
                                IT *)
{
  // Typedefs for image types
  typedef itk::Image< IT, 3 >     ImageType;
  typedef itk::Image< unsigned int, 3 >     OutputImageType;
  typedef itk::Image< unsigned char, 3 >     BinaryImageType;
  // Typedef for importing image data from VolView data
  typedef itk::ImportImageFilter< IT, 3 > ImportFilterType;

  typedef itk::ConnectedComponentImageFilter< ImageType, OutputImageType > ConCompFilterType;
  typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;

  int imgSize = info->InputVolumeDimensions[0]*info->InputVolumeDimensions[1]*
                info->InputVolumeDimensions[2];

  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
  typename ImportFilterType::Pointer maskImportFilter = ImportFilterType::New();

  typename ImportFilterType::IndexType start;
  start.Fill(0);

  typename ImportFilterType::SizeType size;
  size[0] = info->InputVolumeDimensions[0];
  size[1] = info->InputVolumeDimensions[1];
  size[2] = info->InputVolumeDimensions[2];

  typename ImportFilterType::RegionType region;
  region.SetIndex(start);
  region.SetSize(size);

  importFilter->SetRegion(region);

  double spacing[3];
  spacing[0] = info->InputVolumeSpacing[0];
  spacing[1] = info->InputVolumeSpacing[1];
  spacing[2] = info->InputVolumeSpacing[2];
  double origin[3];
  origin[0] = info->InputVolumeOrigin[0];
  origin[1] = info->InputVolumeOrigin[1];
  origin[2] = info->InputVolumeOrigin[2];

  importFilter->SetSpacing(spacing);
  importFilter->SetOrigin(origin);

  // Set data for full image import filter
  importFilter->SetImportPointer((IT*)(pds->inData), imgSize, false);


  // Check for masking data
  if (info->InputVolume2NumberOfComponents > 0) {
    typename ImportFilterType::IndexType startMask;
    startMask.Fill(0);

    typename ImportFilterType::SizeType sizeMask;
    sizeMask[0] = info->InputVolume2Dimensions[0];
    sizeMask[1] = info->InputVolume2Dimensions[1];
    sizeMask[2] = info->InputVolume2Dimensions[2];

    typename ImportFilterType::RegionType regionMask;
    regionMask.SetIndex(startMask);
    regionMask.SetSize(sizeMask);
   
    double spacingMask[3];
    spacingMask[0] = info->InputVolume2Spacing[0];
    spacingMask[1] = info->InputVolume2Spacing[1];
    spacingMask[2] = info->InputVolume2Spacing[2];
    double originMask[3];
    originMask[0] = info->InputVolume2Origin[0];
    originMask[1] = info->InputVolume2Origin[1];
    originMask[2] = info->InputVolume2Origin[2];

    maskImportFilter->SetSpacing(spacingMask);
    maskImportFilter->SetOrigin(originMask);
    maskImportFilter->SetRegion(regionMask);

    int imgSizeMask = info->InputVolume2Dimensions[0]*info->InputVolume2Dimensions[1]*
              info->InputVolume2Dimensions[2];

    // Set data for mask image import filter
    maskImportFilter->SetImportPointer((IT*)(pds->inData2), imgSizeMask, false);

  }

  info->UpdateProgress(info, (float)0.3, "Running filter");

  typename ConCompFilterType::Pointer conCompFilter = ConCompFilterType::New();
  typename RelabelType::Pointer relabelFilter = RelabelType::New();

  conCompFilter->SetInput(importFilter->GetOutput());
  conCompFilter->FullyConnectedOn();
  if (info->InputVolume2NumberOfComponents > 0) {
    conCompFilter->SetMaskImage(maskImportFilter->GetOutput());
  }

  conCompFilter->Update();

  relabelFilter->SetInput(conCompFilter->GetOutput());

  relabelFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
          static_cast< unsigned int * >( pds->outData ),
          imgSize, false);
  relabelFilter->Update();

  std::cout << "There are " << relabelFilter->GetNumberOfObjects() << " objects" <<
    " in this image\n";

  info->UpdateProgress(info,(float)1.0,"Processing Complete");

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
    vtkTemplateMacro3(vvConnectedComponentTemplate, info, pds, 
                      static_cast<VTK_TT *>(0));
  }
  return 0;

}


static int UpdateGUI(void *inf)
{
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  info->OutputVolumeScalarType = VTK_UNSIGNED_INT /*info->InputVolumeScalarType*/;

  info->OutputVolumeNumberOfComponents = 1;

  for (int i = 0; i < 3; i++) {
    info->OutputVolumeDimensions[i] = info->InputVolumeDimensions[i];
    info->OutputVolumeSpacing[i] = info->InputVolumeSpacing[i];
    info->OutputVolumeOrigin[i] = info->InputVolumeOrigin[i];
  }

  return 1;
}


extern "C" {
  
void VV_PLUGIN_EXPORT vvConnectedComponentInit(vtkVVPluginInfo *info)
{
  vvPluginVersionCheck();

  // setup information that never changes
  info->ProcessData = ProcessData;
  info->UpdateGUI   = UpdateGUI;
  info->SetProperty(info, VVP_NAME, "Connected Component");
  info->SetProperty(info, VVP_GROUP, "Utility");
  info->SetProperty(info, VVP_TERSE_DOCUMENTATION,
                                    "Calculates connected components of binarized image.");
  info->SetProperty(info, VVP_FULL_DOCUMENTATION,
    "Simple use of ITK connected component filter.");

  info->SetProperty(info, VVP_SUPPORTS_IN_PLACE_PROCESSING, "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_PIECES,   "0");
  info->SetProperty(info, VVP_NUMBER_OF_GUI_ITEMS,          "0");
  info->SetProperty(info, VVP_REQUIRED_Z_OVERLAP,           "0");
  info->SetProperty(info, VVP_PER_VOXEL_MEMORY_REQUIRED,    "64");
  info->SetProperty(info, VVP_REQUIRES_SECOND_INPUT,        "1");
  info->SetProperty(info, VVP_SECOND_INPUT_OPTIONAL,        "1");
  info->SetProperty(info, VVP_REQUIRES_SERIES_INPUT,        "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_SERIES_BY_VOLUMES, "0");
  info->SetProperty(info, VVP_PRODUCES_OUTPUT_SERIES, "0");
  info->SetProperty(info, VVP_PRODUCES_PLOTTING_OUTPUT, "0");
}

}
