/*=========================================================================

  Copyright (c) Yue Sheng, University of Manchester and David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file vvQuan3D.cpp
\brief This is a VolView plugin to calculate several characteristics from a 
labelled image. It uses the CCPi class CCPiQuantification3D to do all the work.

The following characteristics are calculated:
\li Volume by voxel counts
\li Equivalent sphere diameter by voxel counts
\li Bounding box diagonal
\li Principal Component Analysis
\li Ellipsoid fitting by PCA
\li Equivalent circle diameter by PCA
\li Isosurface by marching cube
\li Surface area
\li Surface volume
\li Equivalent sphere diameter from surface volume
\li Sphercity
\li Normalised surface area to volume ratio (Radius*Sa/Vol)
*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "Quan3D.hpp"
#include "QuanWorker.hpp"


#include "vtkVVPluginAPI.h"

/**
\brief Template function to calculate the data from the image

\c IT represents the input volumes data type (e.g. float, short, etc).

\author David Worth, STFC
\date May 2012

\param info Pointer to object that holds data for this plugin. The voxel size 
can be obtained from it.
\param pds Pointer to object that carries information on data set to be 
processed. Includes actual buffer of voxel data, number of voxels along each
dimension, voxel spacing and voxel type.
*/
template <class IT>
void vvQuan3DTemplate(vtkVVPluginInfo *info,
                      vtkVVProcessDataStruct *pds, 
                      IT *)
{

  CCPiQuantification3D<IT> quan3D;

  info->UpdateProgress(info,0.01,"Initialising...");

  quan3D.Initialise((void*)info, (void*)pds);

  quan3D.CreateVoxelIndexList();

  quan3D.PrepareForQuantification();

  quan3D.PrintSummaryData();  

  quan3D.WriteCSVData("3D_data.csv");

  int totalVoxels = quan3D.GetNumVoxelValues(), n = 0;

  info->UpdateProgress(info,0.05,"Processing...");

  #pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < totalVoxels; i++) {

    CCPiQuantificationWorker<IT> *worker = NULL;
    
    // Do the real work
    #pragma omp critical(nextworker)
    {
      worker = quan3D.GetNextWorker();
    }
    if (worker != NULL) {

      if (0 == worker->Run()) {
        #pragma omp critical(writefile)
        {
          worker->WriteCSVData("3D_data.csv");
        }
      }
      delete worker;
    }
    #pragma omp atomic
    n++;
    #pragma omp critical(progress)
    {
      int count = 40*n/totalVoxels;
      std::cout << '\r';
      for (int j = 1; j < count+1; j++) {
        std::cout << '*';
      }
      for (int j = count+1; j < 41; j++) {
        std::cout << '.';
      }
      std::cout << "|   " << 100*n/totalVoxels << "%";      
    }
  }
  std::cout << std::endl;
/*
  #pragma omp parallel for
  for(int i = 0; i < totalVoxels; i++) {

    #pragma omp flush (abort)
    if (!abort) {

      // See if we should abort
      if (atoi(info->GetProperty(info,VVP_ABORT_PROCESSING)) == 1) {
        abort = true;
        #pragma omp flush (abort)
      }

      #pragma omp critical(UpdateProgress)
      info->UpdateProgress(info,(float)1.0*n/totalVoxels,"Processing...");

      // Do the real work
      CCPiQuantificationWorker<IT> *worker = quan3D.GetWorker(i);
      if (worker != NULL) {
        worker->Run();
        delete worker;
      }

      #pragma omp critical(UpdateProgress)
      ++n;
    }
  }
*/
//  while (quan3D.NextValueQuantification() != 0) {
    /* Show user we're doing something! */
//    info->UpdateProgress(info,(float)1.0*counter/totalVoxels,"Processing...");
    /* See if we should abort */
//    int abort = atoi(info->GetProperty(info,VVP_ABORT_PROCESSING));
//    if (abort) break;
//    counter++;
//  }

  info->UpdateProgress(info,(float)1.0,"Processing Complete");


}

/**
\brief Function called by VolView to run the plug in.

It hands off the actual work to the template function vvQuan3DTemplate.

\param inf Pointer to object that holds data for this plugin. The voxel size 
can be obtained from it.
\param pds Pointer to object that carries information on data set to be 
processed. Includes actual buffer of voxel data, number of voxels along each
dimension, voxel spacing and voxel type.
*/
static int ProcessData(void *inf, vtkVVProcessDataStruct *pds)
{
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  switch (info->InputVolumeScalarType) {
    vtkTemplateMacro3(vvQuan3DTemplate, info, pds, 
                      static_cast<VTK_TT *>(0));
  }
  return 0;
}

/** 
\brief Update the VolView GUI to display user parameters.

Sets one GUI parameter - the physical size of the voxels in the image. It 
gets called prior to the plugin executing.

\param inf Pointer to object that should be modified to set up GUI elements for
this plugin. It also contains details of the input and output images.
*/
static int UpdateGUI(void *inf)
{
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  /* Create required GUI elements here */
  info->SetGUIProperty(info, 0, VVP_GUI_LABEL, "Minimum feature Size");
  info->SetGUIProperty(info, 0, VVP_GUI_TYPE, VVP_GUI_SCALE);
  info->SetGUIProperty(info, 0, VVP_GUI_DEFAULT , "100");
  info->SetGUIProperty(info, 0, VVP_GUI_HELP,
               "Minimum number of voxels for a feature to be analysed");
               
  /* Range for possible output values */
  info->SetGUIProperty(info, 0, VVP_GUI_HINTS , "1 1000 1");


  /* By default the output image's properties match those of the input */
  info->OutputVolumeScalarType = info->InputVolumeScalarType;
  info->OutputVolumeNumberOfComponents = info->InputVolumeNumberOfComponents;
  for (int i = 0; i < 3; i++) {
    info->OutputVolumeDimensions[i] = info->InputVolumeDimensions[i];
    info->OutputVolumeSpacing[i] = info->InputVolumeSpacing[i];
    info->OutputVolumeOrigin[i] = info->InputVolumeOrigin[i];
    }

  return 1;
}

/** 
\brief Initialise this plugin.

Sets the name and group for this plugin in the plugin list shown to the VolView
user and gives some documentation. Also defines properties so VolView can judge
the memory requirements and potential for undoing this plugin.

\param info Pointer to object that should be modified to give details about this
plugin.
*/
extern "C" 
{
  void VV_PLUGIN_EXPORT vvQuan3DInit(vtkVVPluginInfo *info)
  {
    /* Always check the version */
    vvPluginVersionCheck();
    
    /* Set up information that never changes */
    info->ProcessData = ProcessData;
    info->UpdateGUI = UpdateGUI;
    
    /* Set the properties this plugin uses */
    info->SetProperty(info, VVP_NAME, "Quantify 3D");
    info->SetProperty(info, VVP_GROUP, "Quantification");
    
    /* Set terse and full documentation displayed to user */
    info->SetProperty(info, VVP_TERSE_DOCUMENTATION,
     "Quantify several characteristics from a labelled image");
    info->SetProperty(info, VVP_FULL_DOCUMENTATION,
      "Quantify the following characteristics from a labelled image: Volume by \
voxel counts, Equivalent sphere diameter by voxel counts, Bounding box \
diagonal, Principal Component Analysis, Ellipsoid fitting by PCA, Equivalent \
circle diameter by PCA, Isosurface by marching cube, Surface area, Surface \
volume, Equivalent sphere diameter from surface volume, Sphercity, Normalised \
surface area to volume ratio. \nPLEASE NOTE that the Cancel button and progress \
bar do not work for this OpenMP plug-in.");
 
    /* Set these two values to "0" or "1" based on how your plugin
     * handles data all possible combinations of 0 and 1 are valid. */
    info->SetProperty(info, VVP_SUPPORTS_IN_PLACE_PROCESSING, "1");
    info->SetProperty(info, VVP_SUPPORTS_PROCESSING_PIECES,   "0");

    /* Set the number of GUI items used by this plugin */
    info->SetProperty(info, VVP_NUMBER_OF_GUI_ITEMS,          "1");

    info->SetProperty(info, VVP_REQUIRED_Z_OVERLAP,           "0");
    info->SetProperty(info, VVP_REQUIRES_SERIES_INPUT,        "0");
    info->SetProperty(info, VVP_SUPPORTS_PROCESSING_SERIES_BY_VOLUMES, "0");
    info->SetProperty(info, VVP_PRODUCES_OUTPUT_SERIES, "0");
    info->SetProperty(info, VVP_PRODUCES_PLOTTING_OUTPUT, "0");
  }
}



