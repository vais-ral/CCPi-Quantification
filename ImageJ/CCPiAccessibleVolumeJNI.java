/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.0
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class CCPiAccessibleVolumeJNI {
  public final static native long new_MapType__SWIG_0();
  public final static native long new_MapType__SWIG_1(long jarg1, MapType jarg1_);
  public final static native long MapType_size(long jarg1, MapType jarg1_);
  public final static native boolean MapType_empty(long jarg1, MapType jarg1_);
  public final static native void MapType_clear(long jarg1, MapType jarg1_);
  public final static native double MapType_get(long jarg1, MapType jarg1_, double jarg2);
  public final static native void MapType_set(long jarg1, MapType jarg1_, double jarg2, double jarg3);
  public final static native void MapType_del(long jarg1, MapType jarg1_, double jarg2);
  public final static native boolean MapType_has_key(long jarg1, MapType jarg1_, double jarg2);
  public final static native double MapType_getKey(long jarg1, MapType jarg1_, long jarg2);
  public final static native double MapType_getValue(long jarg1, MapType jarg1_, long jarg2);
  public final static native void delete_MapType(long jarg1);
  public final static native long new_CCPiAccessibleVolumeInputImages(int[] jarg1, float[] jarg2, float[] jarg3, long jarg4, CCPiImageDataUnsignedChar jarg4_, long jarg5, CCPiImageDataUnsignedChar jarg5_);
  public final static native void delete_CCPiAccessibleVolumeInputImages(long jarg1);
  public final static native double CCPiAccessibleVolumeInputImages_getScafoldVolume(long jarg1, CCPiAccessibleVolumeInputImages jarg1_);
  public final static native double CCPiAccessibleVolumeInputImages_getScafoldPorosity(long jarg1, CCPiAccessibleVolumeInputImages jarg1_);
  public final static native long[] CCPiAccessibleVolumeInputImages_getDimensions(long jarg1, CCPiAccessibleVolumeInputImages jarg1_);
  public final static native float[] CCPiAccessibleVolumeInputImages_getVoxelSize(long jarg1, CCPiAccessibleVolumeInputImages jarg1_);
  public final static native float[] CCPiAccessibleVolumeInputImages_getOrigin(long jarg1, CCPiAccessibleVolumeInputImages jarg1_);
  public final static native long new_CCPiAccessibleVolumeITKImpl(long jarg1, CCPiAccessibleVolumeInputImages jarg1_, long jarg2, CCPiUserApplicationInterface jarg2_, long jarg3, CCPiImageDataUnsignedChar jarg3_, float jarg4, float jarg5, int jarg6, float jarg7);
  public final static native void delete_CCPiAccessibleVolumeITKImpl(long jarg1);
  public final static native void CCPiAccessibleVolumeITKImpl_Compute(long jarg1, CCPiAccessibleVolumeITKImpl jarg1_);
  public final static native long CCPiAccessibleVolumeITKImpl_GetAccessibleVolume(long jarg1, CCPiAccessibleVolumeITKImpl jarg1_);
  public final static native void CCPiAccessibleVolumeITKImpl_SetOutputImage(long jarg1, CCPiAccessibleVolumeITKImpl jarg1_, long jarg2, CCPiImageDataUnsignedChar jarg2_);
  public final static native long CCPiAccessibleVolumeITKImpl_GetOutputImage(long jarg1, CCPiAccessibleVolumeITKImpl jarg1_);
  public final static native long CCPiAccessibleVolumeITKImpl_GetInputImages(long jarg1, CCPiAccessibleVolumeITKImpl jarg1_);
  public final static native void CCPiUserApplicationInterface_LogMessage(long jarg1, CCPiUserApplicationInterface jarg1_, String jarg2);
  public final static native void CCPiUserApplicationInterface_SetStatusMessage(long jarg1, CCPiUserApplicationInterface jarg1_, String jarg2);
  public final static native void CCPiUserApplicationInterface_SetProgressValue(long jarg1, CCPiUserApplicationInterface jarg1_, float jarg2);
  public final static native boolean CCPiUserApplicationInterface_isCancel(long jarg1, CCPiUserApplicationInterface jarg1_);
  public final static native long new_CCPiUserApplicationInterface();
  public final static native void delete_CCPiUserApplicationInterface(long jarg1);
  public final static native void CCPiUserApplicationInterface_director_connect(CCPiUserApplicationInterface obj, long cptr, boolean mem_own, boolean weak_global);
  public final static native void CCPiUserApplicationInterface_change_ownership(CCPiUserApplicationInterface obj, long cptr, boolean take_or_release);
  public final static native void CCPiConsoleUserInterface_LogMessage(long jarg1, CCPiConsoleUserInterface jarg1_, String jarg2);
  public final static native void CCPiConsoleUserInterface_SetStatusMessage(long jarg1, CCPiConsoleUserInterface jarg1_, String jarg2);
  public final static native void CCPiConsoleUserInterface_SetProgressValue(long jarg1, CCPiConsoleUserInterface jarg1_, float jarg2);
  public final static native boolean CCPiConsoleUserInterface_isCancel(long jarg1, CCPiConsoleUserInterface jarg1_);
  public final static native long new_CCPiConsoleUserInterface();
  public final static native void delete_CCPiConsoleUserInterface(long jarg1);
  public final static native long new_CCPiImageDataUnsignedChar__SWIG_0(short[] jarg1, int[] jarg2, boolean jarg3);
  public final static native long new_CCPiImageDataUnsignedChar__SWIG_1(int[] jarg1);
  public final static native void delete_CCPiImageDataUnsignedChar(long jarg1);
  public final static native short[] CCPiImageDataUnsignedChar_GetImage(long jarg1, CCPiImageDataUnsignedChar jarg1_);
  public final static native int[] CCPiImageDataUnsignedChar_GetDimensions(long jarg1, CCPiImageDataUnsignedChar jarg1_);
  public final static native long CCPiConsoleUserInterface_SWIGUpcast(long jarg1);

  public static void SwigDirector_CCPiUserApplicationInterface_LogMessage(CCPiUserApplicationInterface self, String message) {
    self.LogMessage(message);
  }
  public static void SwigDirector_CCPiUserApplicationInterface_SetStatusMessage(CCPiUserApplicationInterface self, String message) {
    self.SetStatusMessage(message);
  }
  public static void SwigDirector_CCPiUserApplicationInterface_SetProgressValue(CCPiUserApplicationInterface self, float value) {
    self.SetProgressValue(value);
  }
  public static boolean SwigDirector_CCPiUserApplicationInterface_isCancel(CCPiUserApplicationInterface self) {
    return self.isCancel();
  }

  private final static native void swig_module_init();
  static {
    swig_module_init();
  }
}
