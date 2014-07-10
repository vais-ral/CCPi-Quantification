/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class FloatVector {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected FloatVector(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(FloatVector obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_FloatVector(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public FloatVector() {
    this(CCPiAccessibleVolumeJNI.new_FloatVector__SWIG_0(), true);
  }

  public FloatVector(long n) {
    this(CCPiAccessibleVolumeJNI.new_FloatVector__SWIG_1(n), true);
  }

  public long size() {
    return CCPiAccessibleVolumeJNI.FloatVector_size(swigCPtr, this);
  }

  public long capacity() {
    return CCPiAccessibleVolumeJNI.FloatVector_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    CCPiAccessibleVolumeJNI.FloatVector_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return CCPiAccessibleVolumeJNI.FloatVector_isEmpty(swigCPtr, this);
  }

  public void clear() {
    CCPiAccessibleVolumeJNI.FloatVector_clear(swigCPtr, this);
  }

  public void add(float x) {
    CCPiAccessibleVolumeJNI.FloatVector_add(swigCPtr, this, x);
  }

  public float get(int i) {
    return CCPiAccessibleVolumeJNI.FloatVector_get(swigCPtr, this, i);
  }

  public void set(int i, float val) {
    CCPiAccessibleVolumeJNI.FloatVector_set(swigCPtr, this, i, val);
  }

}
