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
#include <set>

#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkHardConnectedComponentImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImportImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkRawImageIO.h"
#include "itkRelabelComponentImageFilter.h"

#include "vtkVVPluginAPI.h"

template <class IT>
void vvAccessibleVolumeTemplate(vtkVVPluginInfo *info,
                                vtkVVProcessDataStruct *pds, 
                                IT *)
{
  // Typedefs for image types
  typedef itk::Image< IT, 3 >     ImageType;
  typedef itk::Image< float, 3 >     DistanceMapImageType;
  typedef itk::Image< unsigned short, 3 >     OutputImageType;
  // Typedef for importing image data from VolView data
  typedef itk::ImportImageFilter< IT, 3 > ImportFilterType;

  // Typedef for distance map filters
  typedef itk::DanielssonDistanceMapImageFilter< ImageType, DistanceMapImageType >   FilterType;
  typedef itk::DanielssonDistanceMapImageFilter< OutputImageType, DistanceMapImageType >   SecondFilterType;
  // Typedef for masking the distance filter output
  typedef itk::MaskImageFilter< DistanceMapImageType, ImageType, DistanceMapImageType > MaskFilterType;

  // Typedefs for calculating the volumes
  typedef itk::BinaryThresholdImageFilter< DistanceMapImageType, OutputImageType > ThresholdFilterType;
  typedef itk::ConnectedComponentImageFilter< OutputImageType, OutputImageType > ConCompFilterType;
  typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
  typedef itk::ImageDuplicator< DistanceMapImageType > DistanceMapDuplicatorType;
  typedef itk::ImageDuplicator< OutputImageType > LabelDuplicatorType;

  // Store results as we go along
  std::ostringstream resultsLog;

  IT *imageData = (IT*)(pds->inData);
  IT *maskData = (IT*)(pds->inData2);
  int imgSize = info->InputVolumeDimensions[0]*info->InputVolumeDimensions[1]*
                info->InputVolumeDimensions[2];
  int sumROI = 0;
  int sumMaskedImg = 0;

  for (int i = 0; i < imgSize; i++) {
    if ( (int)(*maskData) > 0) {
      sumROI++;
      if ( (int)(*imageData) > 0 ) {
        sumMaskedImg++;
      }
    }
    imageData++;
    maskData++;
  }
  // Reset image data pointer
  imageData = (IT*)(pds->inData);
  maskData = (IT*)(pds->inData2);

  std::cout << "Volume of scaffold is " << sumROI-sumMaskedImg << std::endl;
  std::cout << "Scaffold porosity (within ROI) is " << 1.0 - double(sumMaskedImg)/
    double(sumROI) << std::endl;

  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
  typename ImportFilterType::Pointer maskImportFilter = ImportFilterType::New();
  typename FilterType::Pointer distanceMapFilter = FilterType::New();
  //distanceMapFilter->ReleaseDataFlagOn();
  distanceMapFilter->InputIsBinaryOn();

  // Output of import filter is input to distance map filter
  distanceMapFilter->SetInput( importFilter->GetOutput() );

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
  importFilter->SetImportPointer(imageData, imgSize, false);

  maskImportFilter->SetSpacing(spacing);
  maskImportFilter->SetOrigin(origin);
  maskImportFilter->SetRegion(region);

  // Set data for mask image import filter
  maskImportFilter->SetImportPointer(maskData, imgSize, false);

  typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  // Image to be masked is output from distance map
  maskFilter->SetInput1(distanceMapFilter->GetOutput());
  // Mask image is output from mask import filter
  maskFilter->SetInput2(maskImportFilter->GetOutput());

  // Output of mask filter is output image
  typename ImportFilterType::RegionType outputRegion;
  outputRegion.SetIndex(start);
  outputRegion.SetSize(size);
  maskFilter->GetOutput()->SetRegions(outputRegion);
  /* Collect the output of he mask so we can write it to file */
  /*maskFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
          static_cast< unsigned short * >( pds->outData ),
          imgSize, false);
  maskFilter->GetOutput()->Allocate( );*/


  info->UpdateProgress(info, (float)0.05, "Running distance map filter");
  maskFilter->Update();
  std::cout << "Distance map complete\n";

  // Write image to file
  /*itk::RawImageIO<unsigned short,3>::Pointer distMapWriter = itk::RawImageIO<unsigned short,3>::New();
  distMapWriter->SetFileName("DistMap_ROI_512_457.raw");
  for (int i = 0; i < 3; i++) {
    distMapWriter->SetOrigin(i,info->InputVolumeOrigin[i]);
    distMapWriter->SetSpacing(i,info->InputVolumeSpacing[i]);
    distMapWriter->SetDimensions(i,info->InputVolumeDimensions[i]);
  }
  distMapWriter->SetNumberOfComponents(1);
  info->UpdateProgress(info, (float)0.6, "Writing distance map to file");
  distMapWriter->Write(pds->outData);
  std::cout << "File writing complete\n";*/

  // Do the volume calulcations
  float logMin = log(atof(((vtkVVPluginInfo*)info)->GetGUIProperty(info, 0, VVP_GUI_VALUE)));
  float logMax = log(atof(((vtkVVPluginInfo*)info)->GetGUIProperty(info, 1, VVP_GUI_VALUE)));
  int numSpheres = atoi(((vtkVVPluginInfo*)info)->GetGUIProperty(info, 2, VVP_GUI_VALUE));

  float imageResolution = atof(((vtkVVPluginInfo*)info)->GetGUIProperty(info, 3, VVP_GUI_VALUE));

  // From experience first distance map takes 3 times as long as calculation for
  // one sphere diameter so
  float progressIncrement = 1.0/float(numSpheres+3);

  info->UpdateProgress(info, 3.0*progressIncrement, "Start of sphere calculations");

  // Write header for results 
  resultsLog << "Sphere diameter (um), Accessible volume fraction\n";

  //#pragma omp parallel for
  for (unsigned short diamIdx = 0; diamIdx < numSpheres; diamIdx++) {
    float selectedDiameter = exp(logMin + diamIdx*(logMax-logMin)/(numSpheres-1.0));
    float selectedRadius = 0.5*selectedDiameter/imageResolution;
    std::cout << "Selected diameter and radius " << selectedDiameter 
      << " " << selectedRadius << std::endl;

    typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
    typename ConCompFilterType::Pointer conCompFilter = ConCompFilterType::New();
    typename RelabelType::Pointer relabel = RelabelType::New();

    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Calculating threshold");

    // Copy relabelled image out of pipeline for re-use later
    typename DistanceMapDuplicatorType::Pointer distMapDuplicator = DistanceMapDuplicatorType::New();
    distMapDuplicator->SetInputImage(  maskFilter->GetOutput() );
    distMapDuplicator->Update();

    threshold->SetInput(distMapDuplicator->GetOutput());
    threshold->SetInsideValue(itk::NumericTraits<unsigned short>::One);
    threshold->SetOutsideValue(itk::NumericTraits<unsigned short>::Zero);
    threshold->SetLowerThreshold(selectedRadius);
    /*threshold->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    //threshold->Update();
    std::cout << "Threshold calculation complete\n";

    conCompFilter->SetInput(threshold->GetOutput());
    conCompFilter->FullyConnectedOn();
    /*conCompFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    //conCompFilter->Update();
    std::cout << "Connected component calculation complete\n"; 

    relabel->SetInput( conCompFilter->GetOutput() );
    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Calculating components and label map");
    /*relabel->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    relabel->Update();
    std::cout << "Relabelling complete\n";
    std::cout << "There are " << relabel->GetNumberOfObjects() << "objects\n";

    // Copy relabelled image out of pipeline for re-use later
    typename LabelDuplicatorType::Pointer duplicator = LabelDuplicatorType::New();
    duplicator->SetInputImage(  relabel->GetOutput() );
    duplicator->Update();
    
    typedef itk::FlatStructuringElement<3> StructuringElementType;
    typename StructuringElementType::RadiusType  radius;
    radius.Fill(2.5);
    StructuringElementType structuringElement = StructuringElementType::Box(radius);

    typedef itk::GrayscaleDilateImageFilter <OutputImageType, OutputImageType, StructuringElementType>
            GrayscaleDilateImageFilterType;
   
    GrayscaleDilateImageFilterType::Pointer dilateFilter
            = GrayscaleDilateImageFilterType::New();
    dilateFilter->SetInput(relabel->GetOutput());
    dilateFilter->SetKernel(structuringElement);
    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Dilating image");
    /*dilateFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    dilateFilter->Update(); 

    itk::ImageRegionConstIterator< OutputImageType > iter(dilateFilter->GetOutput(),
      dilateFilter->GetOutput()->GetRequestedRegion());
    maskData = (IT*)(pds->inData2);
    std::set<int> selectedGroups;
    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Finding unique labels");
    itk::ImageRegionConstIterator< itk::Image< IT, 3 > > maskIter(maskImportFilter->GetOutput(),
      maskImportFilter->GetOutput()->GetRequestedRegion());
    for (iter.GoToBegin(), maskIter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++maskIter) {
      if ( (iter.Get() != 0) && (maskIter.Get() == 0) )selectedGroups.insert(iter.Get());
      //maskData++;
    }
    for (std::set<int>::iterator setIter = selectedGroups.begin(); setIter != selectedGroups.end(); 
         setIter++){
      std::cout << *setIter << " ";
    }
    std::cout << std::endl;

    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Filtering labelled image");
    itk::ImageRegionIterator< OutputImageType > labelImageIter(duplicator->GetOutput(),
      duplicator->GetOutput()->GetRequestedRegion());
    for (labelImageIter.GoToBegin(); !labelImageIter.IsAtEnd(); ++labelImageIter) {
      if ( selectedGroups.find( labelImageIter.Get() ) != selectedGroups.end() ) {
        labelImageIter.Set(1);
      }
      else {
        labelImageIter.Set(0);
      }
    }

    /*duplicator->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    // Values of duplicator output changed above so better call update again  
    duplicator->Update();

    // Do another distance map
    typename SecondFilterType::Pointer secondDistanceMapFilter = SecondFilterType::New();
    secondDistanceMapFilter->SetInput( duplicator->GetOutput() );
    secondDistanceMapFilter->ReleaseDataFlagOn();
    secondDistanceMapFilter->InputIsBinaryOn(); 

    /*secondDistanceMapFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Calculate distance map of labelled image");

    info->UpdateProgress(info, (3.0+diamIdx)*progressIncrement, "Calculating second distance map and threshold");
    typename ThresholdFilterType::Pointer threshold2 = ThresholdFilterType::New();
    threshold2->SetInput(secondDistanceMapFilter->GetOutput());
    threshold2->SetInsideValue(itk::NumericTraits<unsigned short>::One);
    threshold2->SetOutsideValue(itk::NumericTraits<unsigned short>::Zero);
    threshold2->SetLowerThreshold(0);
    threshold2->SetUpperThreshold(selectedRadius);
    /*threshold2->GetOutput()->GetPixelContainer()->SetImportPointer(
            static_cast< unsigned short * >( pds->outData ),
            imgSize, false);*/
    threshold2->Update();
    unsigned short *outputImage = static_cast< unsigned short * >( pds->outData );
    maskData = (IT*)(pds->inData2);
    itk::ImageRegionConstIterator< OutputImageType > iter2(threshold2->GetOutput(),
      threshold2->GetOutput()->GetRequestedRegion());
    int volPath = 0;
    
    for (iter2.GoToBegin(); !iter2.IsAtEnd(); ++iter2) {
      if (iter2.Get()*((int)(*maskData)) != 0) {
        volPath++;
        *outputImage = diamIdx+1;
      }
      maskData++;
      outputImage++;
    }
    std::cout << "Vol path = " << volPath << std::endl;
    std::cout << "Path volume fraction = " << double(volPath)/double(sumROI-sumMaskedImg) << std::endl;
    resultsLog << selectedDiameter << ", " << double(volPath)/double(sumROI-sumMaskedImg) << std::endl;
  }

  // Create a results file and write to it
  std::ofstream resultsFile("Accessible_Volume.csv", std::ios::trunc);
  resultsFile << resultsLog.str();
  resultsFile.close();

  std::cout << "Processing complete\n";
  info->UpdateProgress(info,(float)1.0,"Processing Complete");

  std::cout.flush();

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
    vtkTemplateMacro3(vvAccessibleVolumeTemplate, info, pds, 
                      static_cast<VTK_TT *>(0));
  }
  return 0;

}


static int UpdateGUI(void *inf)
{
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  info->SetGUIProperty(info, 0, VVP_GUI_LABEL, "Minimum Sphere Diameter (um)");
  info->SetGUIProperty(info, 0, VVP_GUI_TYPE, VVP_GUI_SCALE);
  info->SetGUIProperty(info, 0, VVP_GUI_DEFAULT, "80");
  info->SetGUIProperty(info, 0, VVP_GUI_HELP, "Minimum diameter of sphere used in calculation.");
  info->SetGUIProperty(info, 1, VVP_GUI_LABEL, "Maximum Sphere Diameter (um)");
  info->SetGUIProperty(info, 1, VVP_GUI_TYPE, VVP_GUI_SCALE);
  info->SetGUIProperty(info, 1, VVP_GUI_DEFAULT, "600");
  info->SetGUIProperty(info, 1, VVP_GUI_HELP, "Maximum diameter of sphere used in calculation.");
  info->SetGUIProperty(info, 2, VVP_GUI_LABEL, "Number of sphere diameters to use");
  info->SetGUIProperty(info, 2, VVP_GUI_TYPE, VVP_GUI_SCALE);
  info->SetGUIProperty(info, 2, VVP_GUI_DEFAULT, "11");
  info->SetGUIProperty(info, 2, VVP_GUI_HELP, "Number of sphere diameters to calculate with.");
  info->SetGUIProperty(info, 3, VVP_GUI_LABEL, "Image resolution");
  info->SetGUIProperty(info, 3, VVP_GUI_TYPE, VVP_GUI_SCALE);
  info->SetGUIProperty(info, 3, VVP_GUI_DEFAULT, "9");
  info->SetGUIProperty(info, 3, VVP_GUI_HELP, "Resolution of the image.");
  /* Range for possible output values */
  info->SetGUIProperty(info, 0, VVP_GUI_HINTS , "1 500.0 1.0");
  info->SetGUIProperty(info, 1, VVP_GUI_HINTS , "1 1000.0 1.0");
  info->SetGUIProperty(info, 2, VVP_GUI_HINTS , "1 20 1.0");
  info->SetGUIProperty(info, 3, VVP_GUI_HINTS , "1 1000.0 1.0");

  info->OutputVolumeScalarType = VTK_UNSIGNED_SHORT /*info->InputVolumeScalarType*/;

  info->OutputVolumeNumberOfComponents = 1;

  for (int i = 0; i < 3; i++) {
    info->OutputVolumeDimensions[i] = info->InputVolumeDimensions[i];
    info->OutputVolumeSpacing[i] = info->InputVolumeSpacing[i];
    info->OutputVolumeOrigin[i] = info->InputVolumeOrigin[i];
  }

  return 1;
}


extern "C" {
  
void VV_PLUGIN_EXPORT vvAccessibleVolumeInit(vtkVVPluginInfo *info)
{
  vvPluginVersionCheck();

  // setup information that never changes
  info->ProcessData = ProcessData;
  info->UpdateGUI   = UpdateGUI;
  info->SetProperty(info, VVP_NAME, "Accessible Volume");
  info->SetProperty(info, VVP_GROUP, "Quantification");
  info->SetProperty(info, VVP_TERSE_DOCUMENTATION,
                                    "Calculates accessible volume of binarized image\nwith mask.");
  info->SetProperty(info, VVP_FULL_DOCUMENTATION,
    "Uses the algorithm described in Appendix A of PhD thesis \"Non-destructive quantification of tissue scaffolds and augmentation implants using X-ray microtomography\" by Sheng Yue, Department of Materials, Imperial College London.");

  info->SetProperty(info, VVP_SUPPORTS_IN_PLACE_PROCESSING, "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_PIECES,   "0");
  info->SetProperty(info, VVP_NUMBER_OF_GUI_ITEMS,          "4");
  info->SetProperty(info, VVP_REQUIRED_Z_OVERLAP,           "0");
  info->SetProperty(info, VVP_PER_VOXEL_MEMORY_REQUIRED,    "8");
  info->SetProperty(info, VVP_REQUIRES_SECOND_INPUT,        "1");
  info->SetProperty(info, VVP_SECOND_INPUT_OPTIONAL,        "1");
  info->SetProperty(info, VVP_REQUIRES_SERIES_INPUT,        "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_SERIES_BY_VOLUMES, "0");
  info->SetProperty(info, VVP_PRODUCES_OUTPUT_SERIES, "0");
  info->SetProperty(info, VVP_PRODUCES_PLOTTING_OUTPUT, "0");
}

}
