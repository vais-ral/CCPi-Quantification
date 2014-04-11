// AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY.  Place custom code in custominit.h.
#include "CCPiParticleTracking.h"
#include "CCPiLabelQuantification.h"
#include "CCPiKMeansFilterITK.h"
#include "CCPiComputeThreshold.h"
#include "CCPiAccessibleVolume.h"


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void amirapackage_tutorial_init()
{
    static bool isInitialized = false;
    if (isInitialized)
      return;

    isInitialized = true;

    CCPiParticleTracking::initClass();
    CCPiLabelQuantification::initClass();
    CCPiKMeansFilterITK::initClass();
    CCPiComputeThreshold::initClass();
    CCPiAccessibleVolume::initClass();

}


extern "C"
#ifdef _WIN32
__declspec(dllexport)
#endif
void amirapackage_tutorial_finish()
{
    static bool isFinished = false;
    if (isFinished)
      return;

    isFinished = true;


    CCPiAccessibleVolume::exitClass();
    CCPiComputeThreshold::exitClass();
    CCPiKMeansFilterITK::exitClass();
    CCPiLabelQuantification::exitClass();
    CCPiParticleTracking::exitClass();
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
        amirapackage_tutorial_init();
        break;
    case DLL_PROCESS_DETACH:
        amirapackage_tutorial_finish();
        break;
    default:
        ;
    }
    return true;
}
#endif

#if defined(__GNUC__)
void __attribute__((constructor)) soconstructor_tutorial() {
    amirapackage_tutorial_init();
}

void __attribute__((destructor)) sodestructor_tutorial() {
    amirapackage_tutorial_finish();
}
#endif
