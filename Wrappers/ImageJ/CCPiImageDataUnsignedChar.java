/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class CCPiImageDataUnsignedChar {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected CCPiImageDataUnsignedChar(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(CCPiImageDataUnsignedChar obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_CCPiImageDataUnsignedChar(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public CCPiImageDataUnsignedChar(short[] data, int[] dims, boolean deepCopy) {
    this(CCPiAccessibleVolumeJNI.new_CCPiImageDataUnsignedChar__SWIG_0(data, dims, deepCopy), true);
  }

  public CCPiImageDataUnsignedChar(int[] dims) {
    this(CCPiAccessibleVolumeJNI.new_CCPiImageDataUnsignedChar__SWIG_1(dims), true);
  }

  public short[] GetImage() {
    return CCPiAccessibleVolumeJNI.CCPiImageDataUnsignedChar_GetImage(swigCPtr, this);
  }

  public int[] GetDimensions() {
    return CCPiAccessibleVolumeJNI.CCPiImageDataUnsignedChar_GetDimensions(swigCPtr, this);
  }

}
