/* File : CCPi.i */
%module(directors="1") CCPi

%include "std_map.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_list.i"

%{
#include <map>
#include "../Core/CCPiSimpleHistogramThresholdingITKImpl.h"
%}

namespace std {
	%template(FloatVector) vector<float>;
	%template(StringVector) vector<string>;
	%template(IntList) list<int>;
}

%include "../Core/CCPiSimpleHistogramThresholdingITKImpl.h"

%template(CCPiSimpleHistogramThresholdingITKImplUnsignedChar) CCPiSimpleHistogramThresholdingITKImpl<unsigned char>;

