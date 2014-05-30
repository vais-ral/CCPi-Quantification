/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class MapDoubleDouble implements Iterable<Double> {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected MapDoubleDouble(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(MapDoubleDouble obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_MapDoubleDouble(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public MapDoubleDouble() {
    this(CCPiAccessibleVolumeJNI.new_MapDoubleDouble__SWIG_0(), true);
  }

  public MapDoubleDouble(MapDoubleDouble arg0) {
    this(CCPiAccessibleVolumeJNI.new_MapDoubleDouble__SWIG_1(MapDoubleDouble.getCPtr(arg0), arg0), true);
  }

  public long size() {
    return CCPiAccessibleVolumeJNI.MapDoubleDouble_size(swigCPtr, this);
  }

  public boolean empty() {
    return CCPiAccessibleVolumeJNI.MapDoubleDouble_empty(swigCPtr, this);
  }

  public void clear() {
    CCPiAccessibleVolumeJNI.MapDoubleDouble_clear(swigCPtr, this);
  }

  public double get(double key) {
    return CCPiAccessibleVolumeJNI.MapDoubleDouble_get(swigCPtr, this, key);
  }

  public void set(double key, double x) {
    CCPiAccessibleVolumeJNI.MapDoubleDouble_set(swigCPtr, this, key, x);
  }

  public void del(double key) {
    CCPiAccessibleVolumeJNI.MapDoubleDouble_del(swigCPtr, this, key);
  }

  public boolean has_key(double key) {
    return CCPiAccessibleVolumeJNI.MapDoubleDouble_has_key(swigCPtr, this, key);
  }

  public MapIterator iterator() {
    long cPtr = CCPiAccessibleVolumeJNI.MapDoubleDouble_iterator(swigCPtr, this);
    return (cPtr == 0) ? null : new MapIterator(cPtr, true);
  }

}