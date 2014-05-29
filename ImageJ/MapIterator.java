/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class MapIterator implements java.util.Iterator<Double> {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected MapIterator(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(MapIterator obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_MapIterator(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public void remove() throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }

  public Double next() throws java.util.NoSuchElementException {
    if (!hasNext()) {
      throw new java.util.NoSuchElementException();
    }

    return nextImpl();
  }

  public MapIterator(MapDoubleDouble m) {
    this(CCPiAccessibleVolumeJNI.new_MapIterator(MapDoubleDouble.getCPtr(m), m), true);
  }

  public boolean hasNext() {
    return CCPiAccessibleVolumeJNI.MapIterator_hasNext(swigCPtr, this);
  }

  private double nextImpl() {
    return CCPiAccessibleVolumeJNI.MapIterator_nextImpl(swigCPtr, this);
  }

}
