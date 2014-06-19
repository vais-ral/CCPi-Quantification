/* File : CCPiAccessibleVolume.i */
%module(directors="1") CCPiAccessibleVolume

%include "arrays_java.i"
%include "std_map.i"
%include "std_string.i"
%apply int[] {int *};
%apply unsigned int[] {unsigned int *};
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
#include "..\Core\CCPiAccessibleVolumeInputImages.h"
#include "..\Core\CCPiAccessibleVolumeITKImpl.h"
#include "..\Core\CCPiUserApplicationInterface.h"
#include "..\Core\CCPiConsoleUserInterface.h"
%}

%include "..\Core\CCPiAccessibleVolumeInputImages.h"
%include "..\Core\CCPiAccessibleVolumeITKImpl.h"
%include "..\Core\CCPiUserApplicationInterface.h"
%include "..\Core\CCPiConsoleUserInterface.h"


/*  long size = arg1->GetInputImages()->getDimensions()[0]*arg1->GetInputImages()->getDimensions()[1]*arg1->GetInputImages()->getDimensions()[2];*/