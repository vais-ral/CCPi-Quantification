// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.  Place custom code in custominit.h.
void mcExitClass_CCPiSimpleHistogramThresholding();
void mcInitClass_CCPiSimpleHistogramThresholding();
void mcExitClass_CCPiParticleTracking();
void mcInitClass_CCPiParticleTracking();
void mcExitClass_CCPiLabelQuantification();
void mcInitClass_CCPiLabelQuantification();
void mcExitClass_CCPiKMeansFilterITK();
void mcInitClass_CCPiKMeansFilterITK();
void mcExitClass_CCPiAccessibleVolume();
void mcInitClass_CCPiAccessibleVolume();


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void amirapackage_CCPi_init()
{
    static bool isInitialized = false;
    if (isInitialized)
      return;

    isInitialized = true;

    mcInitClass_CCPiSimpleHistogramThresholding();
    mcInitClass_CCPiParticleTracking();
    mcInitClass_CCPiLabelQuantification();
    mcInitClass_CCPiKMeansFilterITK();
    mcInitClass_CCPiAccessibleVolume();

}


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void amirapackage_CCPi_finish()
{
    static bool isFinished = false;
    if (isFinished)
      return;

    isFinished = true;


    mcExitClass_CCPiAccessibleVolume();
    mcExitClass_CCPiKMeansFilterITK();
    mcExitClass_CCPiLabelQuantification();
    mcExitClass_CCPiParticleTracking();
    mcExitClass_CCPiSimpleHistogramThresholding();
}

#if defined(_WIN32)
#  include <windows.h>
BOOL WINAPI DllMain(
    __in  HINSTANCE hinstDLL,
    __in  DWORD fdwReason,
    __in  LPVOID lpvReserved
    )
{
    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:
        amirapackage_CCPi_init();
        break;
    case DLL_PROCESS_DETACH:
        amirapackage_CCPi_finish();
        break;
    default:
        ;
    }
    return true;
}
#endif

#if defined(__GNUC__)
void __attribute__((constructor)) soconstructor_CCPi() {
    amirapackage_CCPi_init();
}

void __attribute__((destructor)) sodestructor_CCPi() {
    amirapackage_CCPi_finish();
}
#endif
