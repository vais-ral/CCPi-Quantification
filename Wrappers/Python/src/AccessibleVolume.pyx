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
from libcpp.list cimport list
from libcpp.map cimport map

import numpy as np
cimport numpy as np


cdef extern:
	double* AccessibleVolumeProcess(unsigned char* volumeData, unsigned char* maskedVolumeData, long inputDataDims[3], float origin[3], float voxelSize[3], float sphereDiameterRangeMin, float sphereDiameterRangeMax, int numberOfSpheres, float imageResolution)
	

def AccessibleVolume(np.ndarray[np.uint8_t, ndim=3, mode="c"] volumeData, np.ndarray[np.uint8_t, ndim=3, mode="c"] maskedVolumeData, np.ndarray[np.float32_t, ndim=1] nd_origin, np.ndarray[np.float32_t, ndim=1] nd_voxelsize,float sphereDiameterRangeMin, float sphereDiameterRangeMax, int numberOfSpheres, float imageResolution):
	cdef long dims[3]
	dims[0] = volumeData.shape[0]
	dims[1] = volumeData.shape[1]
	dims[2] = volumeData.shape[2]
	cdef float origin[3]
	cdef float voxelsize[3]	
	origin[0] = nd_origin[0]
	origin[1] = nd_origin[1]
	origin[2] = nd_origin[2]
	voxelsize[0] = nd_voxelsize[0]
	voxelsize[1] = nd_voxelsize[1]
	voxelsize[2] = nd_voxelsize[2]
	cdef double *result
	result = AccessibleVolumeProcess(&volumeData[0,0,0], &maskedVolumeData[0,0,0], dims, origin, voxelsize, sphereDiameterRangeMin, sphereDiameterRangeMax, numberOfSpheres, imageResolution)
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] returnValues = np.empty([numberOfSpheres, 2], dtype='float32')
	for i in range(0, numberOfSpheres):
		returnValues[i, 0] = result[i]
		returnValues[i, 1] = result[i+numberOfSpheres]
	return returnValues

		