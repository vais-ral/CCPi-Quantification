/* File : CCPiAccessibleVolume.i */
%module(directors="1") CCPiAccessibleVolume

%include "arrays_java.i"
%include "std_map.i"
%include "std_string.i"
%apply int[] {int *};
%apply unsigned int[] {unsigned int *};
%apply float[] {float *};
%apply unsigned char[] {unsigned char *};

%typemap(javainterfaces) MapIterator "java.util.Iterator<Double>"
%typemap(javacode) MapIterator %{
  public void remove() throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }

  public Double next() throws java.util.NoSuchElementException {
    if (!hasNext()) {
      throw new java.util.NoSuchElementException();
    }

    return nextImpl();
  }
%}
	
%javamethodmodifiers MapIterator::nextImpl "private";
%inline %{
  struct MapIterator {
    typedef std::map<double,double> map_t;
    MapIterator(const map_t& m) : it(m.begin()), map(m) {}
    bool hasNext() const {
      return it != map.end();
    }

    const double& nextImpl() {
      const std::pair<double,double>& ret = *it++;
      return ret.second;
    }
  private:
    map_t::const_iterator it;
    const map_t& map;    
  };
%}
	
%typemap(javainterfaces) std::map<double,double> "Iterable<Double>"

%newobject std::map<double,double>::iterator() const;
%extend std::map<double,double> {
  MapIterator *iterator() const {
    return new MapIterator(*$self);
  }
}

%template (MapDoubleDouble) std::map<double,double>;

/*
%typemap(javain) CCPiAccessibleVolumeInputImages *input "getCPtrAndAddReferenceInputImages($javainput)"
%typemap(javain) CCPiUserApplicationInterface *userAppInterface "getCPtrAndAddReferenceUI($javainput)"
%typemap(javacode) CCPiAccessibleVolumeITKImpl %{
  // Ensure that the GC doesn't collect any element set from Java
  // as the underlying C++ class stores a shallow copy
  private CCPiAccessibleVolumeInputImages inputDataReference;
  private CCPiUserApplicationInterface uiReference;
  private long getCPtrAndAddReferenceInputImages(CCPiAccessibleVolumeInputImages inputData) {
    inputDataReference = inputData;
    return CCPiAccessibleVolumeInputImages.getCPtr(inputData);
  }
  private long getCPtrAndAddReferenceUI(CCPiUserApplicationInterface ui) {
    uiReference = ui;
    return CCPiUserApplicationInterface.getCPtr(ui);
  }
%}
*/

%feature("director") CCPiUserApplicationInterface;

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