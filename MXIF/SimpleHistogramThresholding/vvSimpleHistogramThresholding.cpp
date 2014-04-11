/*=========================================================================

=========================================================================*/
/** \file vvSimpleHistogramThresholding.cpp
 *
 * This plugin performs a simple thresholding of an image. It is designed for
 * images which have histogram with two distinct peaks and sets the threshold
 * value to be the average of the two peak intensity values.
 * It will not work well for other images!
 */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <set>

#include "itkScalarImageToHistogramGenerator.h"
#include "itkImage.h"
#include "itkImportImageFilter.h"

#include "vtkVVPluginAPI.h"

template <class IT>
void vvSimpleHistogramThresholdingTemplate(vtkVVPluginInfo *info,
                               vtkVVProcessDataStruct *pds, 
                               IT *)
{
  // Typedefs for image types
  typedef itk::Image< IT, 3 >     ImageType;
  typedef itk::Image< unsigned char, 3 >     OutputImageType;
  // Typedef for importing image data from VolView data
  typedef itk::ImportImageFilter< IT, 3 > ImportFilterType;

  // Typedef for filter
  typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > HistogramGeneratorType;


  IT *imageData = (IT*)(pds->inData);
  int imgSize = info->InputVolumeDimensions[0]*info->InputVolumeDimensions[1]*
                info->InputVolumeDimensions[2];

  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
  typename HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();


  // Output of import filter is input to distance map filter
  histogramGenerator->SetInput( importFilter->GetOutput() );

  typename ImportFilterType::IndexType start;
  start.Fill(0);

  typename ImportFilterType::SizeType size;
  size[0] = info->InputVolumeDimensions[0];
  size[1] = info->InputVolumeDimensions[1];
  size[2] = info->InputVolumeDimensions[2];

  typename ImportFilterType::RegionType region;
  region.SetIndex(start);
  region.SetSize(size);
 
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
  importFilter->SetRegion(region);

  // Set data for full image import filter
  importFilter->SetImportPointer(imageData+size[0] * size[1]*pds->StartSlice, imgSize, false);
  importFilter->Update();

  info->UpdateProgress(info, (float)0.3, "Calculating Histogram");
  // Number of bins = range of input image so 1 bin covers 1 intensity
  const unsigned int numBins = info->InputVolumeScalarRange[1]-info->InputVolumeScalarRange[0];
  histogramGenerator->SetNumberOfBins(numBins);
  histogramGenerator->SetMarginalScale( 10.0 );

  histogramGenerator->SetHistogramMin(info->InputVolumeScalarRange[0]-0.5);
  histogramGenerator->SetHistogramMax(info->InputVolumeScalarRange[1]+0.5);
  histogramGenerator->Compute();  
  
  typedef typename HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = histogramGenerator->GetOutput();

  const unsigned int histogramSize = histogram->Size();

  std::cout << "Histogram size " << histogramSize << std::endl;

  // Want to find first and last peak so have array size 2 to store bin numbers
  bool goingUp = false;
  int peaks[2];
  int peakId = 0;
  
  for( unsigned int bin=0; bin < histogramSize-1; bin++ ) {
    std::cout << "bin = " << bin << " frequency = ";
    std::cout << histogram->GetFrequency( bin, 0 ) << std::endl;
    if ( histogram->GetFrequency( bin+1, 0 ) < histogram->GetFrequency( bin, 0 ) ) {
      if (goingUp == true ) {
        peaks[peakId] = bin;
        peakId = 1;
        goingUp = false;
        std::cout << "Peak at bin " << bin << std::endl;
      }
    }
    else {
      goingUp = true;
    }
  }
  
  std::cout << "First and last peaks found at bin " << peaks[0] << " and " << peaks[1] << std::endl;
  std::cout << "Correpsonding intensity values are " << info->InputVolumeScalarRange[0] + peaks[0] <<
               " and " << info->InputVolumeScalarRange[0] + peaks[1] << std::endl;
  
  
  /*kmeansFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned char * >( pds->outData ),
            imgSize, false);*/




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
    vtkTemplateMacro3(vvSimpleHistogramThresholdingTemplate, info, pds, 
                      static_cast<VTK_TT *>(0));
  }
  return 0;

}


static int UpdateGUI(void *inf)
{
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  info->OutputVolumeScalarType = VTK_UNSIGNED_CHAR;

  info->OutputVolumeNumberOfComponents = 1;

  for (int i = 0; i < 3; i++) {
    info->OutputVolumeDimensions[i] = info->InputVolumeDimensions[i];
    info->OutputVolumeSpacing[i] = info->InputVolumeSpacing[i];
    info->OutputVolumeOrigin[i] = info->InputVolumeOrigin[i];
  }

  return 1;
}


extern "C" {
  
void VV_PLUGIN_EXPORT vvSimpleHistogramThresholdingInit(vtkVVPluginInfo *info)
{
  vvPluginVersionCheck();

  // setup information that never changes
  info->ProcessData = ProcessData;
  info->UpdateGUI   = UpdateGUI;
  info->SetProperty(info, VVP_NAME, "Simple Histogram Thresholding");
  info->SetProperty(info, VVP_GROUP, "Segmentation - Statistics");
  info->SetProperty(info, VVP_TERSE_DOCUMENTATION,
                                    "Thresholds image based on intensity histogram.");
  info->SetProperty(info, VVP_FULL_DOCUMENTATION,
    "It is designed for images which have histogram with two distinct peaks \
    and sets the threshold value to be the average of the two peak intensity values. \
    It will not work well for other images!.");

  info->SetProperty(info, VVP_SUPPORTS_IN_PLACE_PROCESSING, "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_PIECES,   "0");
  info->SetProperty(info, VVP_NUMBER_OF_GUI_ITEMS,          "0");
  info->SetProperty(info, VVP_REQUIRED_Z_OVERLAP,           "0");
  info->SetProperty(info, VVP_PER_VOXEL_MEMORY_REQUIRED,    "2");
  info->SetProperty(info, VVP_REQUIRES_SERIES_INPUT,        "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_SERIES_BY_VOLUMES, "0");
  info->SetProperty(info, VVP_PRODUCES_OUTPUT_SERIES, "0");
  info->SetProperty(info, VVP_PRODUCES_PLOTTING_OUTPUT, "0");
}

}
