#ifndef CCPIDEFINES_H
#define CCPIDEFINES_H

#if defined(_WIN32) || defined(__WIN32__)
  #if defined(CCPiCore_EXPORTS) // add by CMake 
    #define  CCPI_EXPORT extern "C" __declspec(dllexport)
  #else
    #define  CCPI_EXPORT extern "C" __declspec(dllimport)
  #endif /* CCPi_EXPORTS */
#elif defined(linux) || defined(__linux)
 #define CCPI_EXPORT
#endif

#endif