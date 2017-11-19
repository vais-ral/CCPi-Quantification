/**
Copyright 2017 CCPi

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

AccessibleVolume_c.cpp
Accessible Volume algorithm
Author: Mr. Srikanth Nagella 
*/

#ifdef WIN32
#include <iostream>
#include <boost/shared_ptr.hpp>
#endif

#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiAccessibleVolumeITKImpl.h"
#include "CCPiConsoleUserInterface.h"

#include "vtkType.h"

extern "C" {
double* AccessibleVolumeProcess(unsigned char* volumeData, unsigned char* maskedVolumeData, long inputDataDims[3], float origin[3], float voxelSize[3], float sphereDiameterRangeMin, float sphereDiameterRangeMax, int numberOfSpheres, float imageResolution)
{	
	CCPiImageDataUnsignedChar* inputData = new CCPiImageDataUnsignedChar(volumeData, inputDataDims, true);
	CCPiImageDataUnsignedChar* inputMaskData = new CCPiImageDataUnsignedChar(maskedVolumeData, inputDataDims, true);
	int dims[3];
	dims[0]=inputDataDims[0];dims[1]=inputDataDims[1];dims[2]=inputDataDims[2];
	CCPiAccessibleVolumeInputImages input(dims, voxelSize, origin, inputData, inputMaskData);
	CCPiConsoleUserInterface ui;
	CCPiAccessibleVolumeITKImpl avii(&input, &ui, NULL, sphereDiameterRangeMin, sphereDiameterRangeMax, numberOfSpheres, imageResolution);
	avii.Compute();
	std::map<double,double> outputAV=avii.GetAccessibleVolume();
	double* result=new double[numberOfSpheres*2];
	int i=0;
	for(std::map<double,double>::iterator it = outputAV.begin(); it != outputAV.end(); ++it, i++) {
		result[i]=it->first;
		result[i+numberOfSpheres]=it->second;
	}	
	return result;
}
}