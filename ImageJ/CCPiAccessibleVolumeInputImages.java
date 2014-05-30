/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class CCPiAccessibleVolumeInputImages {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected CCPiAccessibleVolumeInputImages(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(CCPiAccessibleVolumeInputImages obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CCPiAccessibleVolumeJNI.delete_CCPiAccessibleVolumeInputImages(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public CCPiAccessibleVolumeInputImages(int[] volumeDims, float[] voxelSize, float[] origin, short[] volumeData, short[] maskedVolumeData) {
    this(CCPiAccessibleVolumeJNI.new_CCPiAccessibleVolumeInputImages(volumeDims, voxelSize, origin, volumeData, maskedVolumeData), true);
  }

  public double getScafoldVolume() {
    return CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_getScafoldVolume(swigCPtr, this);
  }

  public double getScafoldPorosity() {
    return CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_getScafoldPorosity(swigCPtr, this);
  }

  public long[] getDimensions() {
    return CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_getDimensions(swigCPtr, this);
  }

  public float[] getVoxelSize() {
    return CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_getVoxelSize(swigCPtr, this);
  }

  public float[] getOrigin() {
    return CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_getOrigin(swigCPtr, this);
  }

  public SWIGTYPE_p_itk__ImageT_unsigned_char_3_t__Pointer GetVolumeMaskData() {
    return new SWIGTYPE_p_itk__ImageT_unsigned_char_3_t__Pointer(CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_GetVolumeMaskData(swigCPtr, this), true);
  }

  public SWIGTYPE_p_itk__ImageT_unsigned_char_3_t__Pointer GetVolumeData() {
    return new SWIGTYPE_p_itk__ImageT_unsigned_char_3_t__Pointer(CCPiAccessibleVolumeJNI.CCPiAccessibleVolumeInputImages_GetVolumeData(swigCPtr, this), true);
  }

}