/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class MapType {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected MapType(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(MapType obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_MapType(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public MapType() {
    this(CCPiAccessibleVolumeJNI.new_MapType__SWIG_0(), true);
  }

  public MapType(MapType arg0) {
    this(CCPiAccessibleVolumeJNI.new_MapType__SWIG_1(MapType.getCPtr(arg0), arg0), true);
  }

  public long size() {
    return CCPiAccessibleVolumeJNI.MapType_size(swigCPtr, this);
  }

  public boolean empty() {
    return CCPiAccessibleVolumeJNI.MapType_empty(swigCPtr, this);
  }

  public void clear() {
    CCPiAccessibleVolumeJNI.MapType_clear(swigCPtr, this);
  }

  public double get(double key) {
    return CCPiAccessibleVolumeJNI.MapType_get(swigCPtr, this, key);
  }

  public void set(double key, double x) {
    CCPiAccessibleVolumeJNI.MapType_set(swigCPtr, this, key, x);
  }

  public void del(double key) {
    CCPiAccessibleVolumeJNI.MapType_del(swigCPtr, this, key);
  }

  public boolean has_key(double key) {
    return CCPiAccessibleVolumeJNI.MapType_has_key(swigCPtr, this, key);
  }

  public double getKey(long index) {
    return CCPiAccessibleVolumeJNI.MapType_getKey(swigCPtr, this, index);
  }

  public double getValue(long index) {
    return CCPiAccessibleVolumeJNI.MapType_getValue(swigCPtr, this, index);
  }

}