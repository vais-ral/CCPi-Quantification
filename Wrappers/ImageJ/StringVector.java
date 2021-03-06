/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class StringVector {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected StringVector(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(StringVector obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_StringVector(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public StringVector() {
    this(CCPiAccessibleVolumeJNI.new_StringVector__SWIG_0(), true);
  }

  public StringVector(long n) {
    this(CCPiAccessibleVolumeJNI.new_StringVector__SWIG_1(n), true);
  }

  public long size() {
    return CCPiAccessibleVolumeJNI.StringVector_size(swigCPtr, this);
  }

  public long capacity() {
    return CCPiAccessibleVolumeJNI.StringVector_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    CCPiAccessibleVolumeJNI.StringVector_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return CCPiAccessibleVolumeJNI.StringVector_isEmpty(swigCPtr, this);
  }

  public void clear() {
    CCPiAccessibleVolumeJNI.StringVector_clear(swigCPtr, this);
  }

  public void add(String x) {
    CCPiAccessibleVolumeJNI.StringVector_add(swigCPtr, this, x);
  }

  public String get(int i) {
    return CCPiAccessibleVolumeJNI.StringVector_get(swigCPtr, this, i);
  }

  public void set(int i, String val) {
    CCPiAccessibleVolumeJNI.StringVector_set(swigCPtr, this, i, val);
  }

}
