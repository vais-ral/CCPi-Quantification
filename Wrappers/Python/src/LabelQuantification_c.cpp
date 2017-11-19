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

label.pyx
Label Quantification algorithm
Author: Mr. Srikanth Nagella 
*/

#ifdef WIN32
#include <iostream>
#include <boost/shared_ptr.hpp>
#endif

#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiAccessibleVolumeITKImpl.h"
#include "CCPiConsoleUserInterface.h"

#include "CCPiLabelQuantificationITKImpl.h"
#include "vtkType.h"

template <class T>
int getVtkType()
{
	if(typeid(T) == typeid(unsigned char))
		return VTK_SIGNED_CHAR;
	else if(typeid(T) == typeid(char))
		return VTK_CHAR;
	else if(typeid(T) == typeid(short))
		return VTK_SHORT;
	else if(typeid(T) == typeid(unsigned short))
		return VTK_UNSIGNED_SHORT;
	else if(typeid(T) == typeid(int))
		return VTK_INT;
	else if(typeid(T) == typeid(unsigned int))
		return VTK_UNSIGNED_INT;
	else if(typeid(T) == typeid(long))
		return VTK_LONG;
	else if(typeid(T) == typeid(unsigned long))
		return VTK_UNSIGNED_LONG;
	else if(typeid(T) == typeid(float))
		return VTK_FLOAT;
	else if(typeid(T) == typeid(double))
		return VTK_DOUBLE;
	return VTK_ID_TYPE;
}

template<class T>
void LabelQuantificationProcess(T* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{	
	CCPiImageData<T> imageData(data, dims, true);
	CCPiConsoleUserInterface ui;	
	CCPiLabelQuantificationITKImpl<T> lq(&imageData, &ui, origin, dims, voxelSize, min, max, minFeatureSize, getVtkType<T>());
	lq.Compute();
	CCPiLabelQuantificationResult* output=lq.GetOutput();
	std::vector<std::string> columnNames = output->GetQuantityNames();
	for(std::vector<std::string>::iterator i= columnNames.begin();i!=columnNames.end();++i)
	{
		featureNames.push_back(*i);
	}
	std::list<int> labelIndexes = output->GetLabelIndexes();
	int xdim = labelIndexes.size();
	int ydim = columnNames.size();
	*nvals = labelIndexes.size();
	*featureValues = new double[labelIndexes.size()*columnNames.size()];
	int i=0;
	for(std::list<int>::iterator row_itr=labelIndexes.begin();row_itr!=labelIndexes.end();row_itr++, i++)
	{
		int j=0;
		for(std::vector<std::string>::iterator column_itr = columnNames.begin();column_itr!=columnNames.end();column_itr++, j++)
		{
			double value = output->GetValue(*column_itr, *row_itr);
			(*featureValues)[i*ydim+j] = value;
		}
	}	
}

extern "C" {
void LabelQuantificationProcessInt(int* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{
	return LabelQuantificationProcess<int>(data, dims, origin, voxelSize, min,max,minFeatureSize, featureNames, featureValues, nvals) ;
}
void LabelQuantificationProcessShort(short* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{
	return LabelQuantificationProcess<short>(data, dims, origin, voxelSize, min,max,minFeatureSize,  featureNames, featureValues, nvals) ;
}
void LabelQuantificationProcessChar(char* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{
	return LabelQuantificationProcess<char>(data, dims, origin, voxelSize, min,max,minFeatureSize, featureNames, featureValues, nvals) ;
}
void LabelQuantificationProcessUInt(unsigned int* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{
	return LabelQuantificationProcess<unsigned int>(data, dims, origin, voxelSize, min,max,minFeatureSize, featureNames, featureValues, nvals) ;
}
void LabelQuantificationProcessUChar(unsigned char* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{
	return LabelQuantificationProcess<unsigned char>(data, dims, origin, voxelSize, min,max,minFeatureSize, featureNames, featureValues, nvals) ;
}
void LabelQuantificationProcessUShort(unsigned short* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, std::vector<std::string>& featureNames, double** featureValues, int* nvals)
{
	return LabelQuantificationProcess<unsigned short>(data, dims, origin, voxelSize, min,max,minFeatureSize, featureNames, featureValues, nvals) ;
}

}