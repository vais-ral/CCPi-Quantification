# distutils: language=c++
"""
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
"""

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
#from libcpp.list cimport list
from libcpp.map cimport map

import numpy as np
cimport numpy as np


cdef extern:
	void LabelQuantificationProcessInt(int* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, vector[string]& featureNames, double** featureValues, int* nvals)
	void LabelQuantificationProcessShort(short* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, vector[string]& featureNames, double** featureValues, int* nvals)
	void LabelQuantificationProcessChar(char* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, vector[string]& featureNames, double** featureValues, int* nvals)
	void LabelQuantificationProcessUInt(unsigned int* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, vector[string]& featureNames, double** featureValues, int* nvals)
	void LabelQuantificationProcessUShort(unsigned short* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, vector[string]& featureNames, double** featureValues, int* nvals)
	void LabelQuantificationProcessUChar(unsigned char* data, long dims[3], float origin[3], float voxelSize[3],float min, float max, float minFeatureSize, vector[string]& featureNames, double** featureValues, int* nvals)
	

def LabelQuantification(np.ndarray volume, np.ndarray nd_origin, np.ndarray nd_voxelsize,float min, float max, float minFeatureSize ):

	# declare a numpy array of raw bytes (unsigned 8-bit integers)
	# and assign it to a view of the input data.
	cdef np.uint8_t[:, :, :] buffer
	buffer = volume.view(np.uint8)
	 
	cdef long dims[3]
	dims[0] = volume.shape[0]
	dims[1] = volume.shape[1]
	dims[2] = volume.shape[2]
	cdef float origin[3]
	cdef float voxelsize[3]	
	origin[0] = nd_origin[0]
	origin[1] = nd_origin[1]
	origin[2] = nd_origin[2]
	voxelsize[0] = nd_voxelsize[0]
	voxelsize[1] = nd_voxelsize[1]
	voxelsize[2] = nd_voxelsize[2]
	cdef vector[string] featureNames
	cdef double* featureValues
	cdef int nvals
	
	# choose the appropriate routine based
	if volume.dtype == np.int32:
		LabelQuantificationProcessInt(<int *>&buffer[0, 0, 0], dims, origin, voxelsize, min, max, minFeatureSize, featureNames, &featureValues, &nvals)
	elif volume.dtype == np.int16:
		LabelQuantificationProcessShort(<short *>&buffer[0, 0, 0], dims, origin, voxelsize, min, max, minFeatureSize, featureNames, &featureValues, &nvals)
	elif volume.dtype == np.int8:
		LabelQuantificationProcessChar(<char *>&buffer[0, 0, 0], dims, origin, voxelsize, min, max, minFeatureSize, featureNames, &featureValues, &nvals)
	elif volume.dtype == np.uint32:
		LabelQuantificationProcessUInt(<unsigned int *>&buffer[0, 0, 0], dims, origin, voxelsize, min, max, minFeatureSize, featureNames, &featureValues, &nvals)
	elif volume.dtype == np.uint16:
		LabelQuantificationProcessUShort(<unsigned short *>&buffer[0, 0, 0], dims, origin, voxelsize, min, max, minFeatureSize, featureNames, &featureValues, &nvals)
	elif volume.dtype == np.uint8:
		LabelQuantificationProcessUChar(<unsigned char *>&buffer[0, 0, 0], dims, origin, voxelsize, min, max, minFeatureSize, featureNames, &featureValues, &nvals)
	else:
		raise ValueError("dtype {0} not supported".format(volume.dtype))
	
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] returnFeatureValues = np.empty([nvals, featureNames.size()], dtype='float32')
	for i in range(0, nvals):
		for j in range(0, featureNames.size()):
			returnFeatureValues[i, j] = featureValues[i*featureNames.size()+j]
	return featureNames, returnFeatureValues

		
