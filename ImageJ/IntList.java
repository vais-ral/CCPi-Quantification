/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class IntList {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected IntList(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(IntList obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_IntList(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public IntList() {
    this(CCPiAccessibleVolumeJNI.new_IntList(), true);
  }

  public long size() {
    return CCPiAccessibleVolumeJNI.IntList_size(swigCPtr, this);
  }

  public boolean isEmpty() {
    return CCPiAccessibleVolumeJNI.IntList_isEmpty(swigCPtr, this);
  }

  public void clear() {
    CCPiAccessibleVolumeJNI.IntList_clear(swigCPtr, this);
  }

  public void add(int x) {
    CCPiAccessibleVolumeJNI.IntList_add(swigCPtr, this, x);
  }

  public int get(int i) {
    return CCPiAccessibleVolumeJNI.IntList_get(swigCPtr, this, i);
  }

}