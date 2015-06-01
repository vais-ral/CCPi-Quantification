/* File : CCPiCore.i */
%module(directors="1") CCPiCore

%include "numpy.i"


%apply unsigned char[] {unsigned char *};


%{
#include <map>
#include "..\Core\CCPiSimpleHistogramThresholdingITKImpl.h"
%}


%include "..\Core\CCPiSimpleHistogramThresholdingITKImpl.h"

%template(CCPiSimpleHistogramThresholdingITKImplUnsignedChar) CCPiSimpleHistogramThresholdingITKImpl<unsigned char>;

