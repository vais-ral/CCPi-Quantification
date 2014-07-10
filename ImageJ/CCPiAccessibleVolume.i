/* File : CCPiAccessibleVolume.i */
%module(directors="1") CCPiAccessibleVolume

%include "arrays_java.i"
%include "std_map.i"
%include "std_string.i"
%include "std_vector.i"

%apply int[] {long *};
%apply int[] {int *};
%apply unsigned int[] {unsigned int * getDimensions()};
%apply float[] {float *getVoxelSize()};
%apply float[] {float *getOrigin()};
%apply float[] {float *};

%apply unsigned char[] {unsigned char *};

%typemap(jstype) std::map<double, double> "java.util.Map<Double,Double>"
%typemap(javaout) std::map<double, double> {
	java.util.Map<Double,Double> result = new java.util.LinkedHashMap<Double,Double>();
	MapType map = new MapType($jnicall,true);
	for(int idx=0;idx<map.size();idx++)
		result.put(map.getKey(idx),map.getValue(idx));
	return result;
}
%template(MapType) std::map<double, double>;

%feature("director") CCPiUserApplicationInterface;

%ignore CCPiAccessibleVolumeInputImages::GetVolumeData;
%ignore CCPiAccessibleVolumeInputImages::GetVolumeMaskData;

%{
#include <map>
#include "..\Core\CCPiImageData.h"
#include "..\Core\CCPiAccessibleVolumeInputImages.h"
#include "..\Core\CCPiAccessibleVolumeITKImpl.h"
#include "..\Core\CCPiUserApplicationInterface.h"
#include "..\Core\CCPiConsoleUserInterface.h"
#include "..\Core\CCPiSimpleHistogramThresholdingITKImpl.h"
%}


%typemap(out) long* GetDimensions {
 $result = SWIG_JavaArrayOutInt(jenv, (int*)$1, 3 /*FillMeInAsSizeCannotBeDeterminedAutomatically*/); 
} 
%typemap(out) unsigned int* getDimensions {
$result = SWIG_JavaArrayOutUint(jenv, (unsigned int*)$1, 3 /*FillMeInAsSizeCannotBeDeterminedAutomatically*/);
}
%typemap(out) float* getVoxelSize {
$result = SWIG_JavaArrayOutFloat(jenv, $1, 3 /*FillMeInAsSizeCannotBeDeterminedAutomatically*/); 
} 
%typemap(out) float* getOrigin {
$result = SWIG_JavaArrayOutFloat(jenv, $1, 3 /*FillMeInAsSizeCannotBeDeterminedAutomatically*/); 
} 

%typemap(out) unsigned char* GetImage {
  long size = arg1->GetDimensions()[0]*arg1->GetDimensions()[1]*arg1->GetDimensions()[2];
  jresult =	SWIG_JavaArrayOutUchar(jenv, (unsigned char *)result, size);
}

namespace std {
	%template(FloatVector) vector<float>;
}

%include "..\Core\CCPiImageData.h"
%include "..\Core\CCPiAccessibleVolumeInputImages.h"
%include "..\Core\CCPiAccessibleVolumeITKImpl.h"
%include "..\Core\CCPiUserApplicationInterface.h"
%include "..\Core\CCPiConsoleUserInterface.h"
%include "..\Core\CCPiSimpleHistogramThresholdingITKImpl.h"

%template(CCPiImageDataUnsignedChar) CCPiImageData<unsigned char>;
%template(CCPiSimpleHistogramThresholdingITKImplUnsignedChar) CCPiSimpleHistogramThresholdingITKImpl<unsigned char>;
