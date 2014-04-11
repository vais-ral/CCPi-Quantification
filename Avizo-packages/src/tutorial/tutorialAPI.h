/***************************************************************
 *
 * DLL export/import definitions for Windows
 *
 ***************************************************************/

#ifndef TUTORIAL_API_H
#define TUTORIAL_API_H

#ifdef HX_OS_WIN
    #ifdef TUTORIAL_EXPORTS
        #define TUTORIAL_API __declspec(dllexport)
    #else
        #define TUTORIAL_API __declspec(dllimport)
    #endif
#else
    #define TUTORIAL_API
#endif

#endif // TUTORIAL_API_H
