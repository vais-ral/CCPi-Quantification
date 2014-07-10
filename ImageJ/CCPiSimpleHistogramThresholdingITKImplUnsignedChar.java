/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class CCPiSimpleHistogramThresholdingITKImplUnsignedChar {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected CCPiSimpleHistogramThresholdingITKImplUnsignedChar(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(CCPiSimpleHistogramThresholdingITKImplUnsignedChar obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_CCPiSimpleHistogramThresholdingITKImplUnsignedChar(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public CCPiSimpleHistogramThresholdingITKImplUnsignedChar(CCPiImageDataUnsignedChar inputImage, int[] volumeDims, float[] voxelSize, float[] origin, float minIntensity, float maxIntensity) {
    this(CCPiAccessibleVolumeJNI.new_CCPiSimpleHistogramThresholdingITKImplUnsignedChar(CCPiImageDataUnsignedChar.getCPtr(inputImage), inputImage, volumeDims, voxelSize, origin, minIntensity, maxIntensity), true);
  }

  public void Compute() {
    CCPiAccessibleVolumeJNI.CCPiSimpleHistogramThresholdingITKImplUnsignedChar_Compute(swigCPtr, this);
  }

  public CCPiImageDataUnsignedChar GetOutputImage() {
    long cPtr = CCPiAccessibleVolumeJNI.CCPiSimpleHistogramThresholdingITKImplUnsignedChar_GetOutputImage__SWIG_0(swigCPtr, this);
    return (cPtr == 0) ? null : new CCPiImageDataUnsignedChar(cPtr, false);
  }

  public void GetOutputImage(CCPiImageDataUnsignedChar returnImage) {
    CCPiAccessibleVolumeJNI.CCPiSimpleHistogramThresholdingITKImplUnsignedChar_GetOutputImage__SWIG_1(swigCPtr, this, CCPiImageDataUnsignedChar.getCPtr(returnImage), returnImage);
  }

  public FloatVector GetPeaks() {
    return new FloatVector(CCPiAccessibleVolumeJNI.CCPiSimpleHistogramThresholdingITKImplUnsignedChar_GetPeaks(swigCPtr, this), true);
  }

}
