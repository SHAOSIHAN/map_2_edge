// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the BASIC_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// BASIC_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifndef _BASIC_COMMON_H_
#define _BASIC_COMMON_H_



#ifdef WIN32
// Disable a warning message about dll
// this is a temporary solution
// http://support.microsoft.com/default.aspx?scid=kb;EN-US;168958
#   pragma warning( disable : 4251 )
#endif

// Win 32 DLL export macros
#ifdef WIN32
# ifdef BASIC_EXPORTS
#   define BASIC_API  __declspec(dllexport)
# else
#   define BASIC_API  __declspec(dllimport)
# endif
#endif // WIN32



extern BASIC_API int nBasic;

BASIC_API int fnBasic(void);



#endif