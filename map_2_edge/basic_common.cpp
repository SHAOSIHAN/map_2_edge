// Basic.cpp : Defines the exported functions for the DLL application.
//

#include "basic_common.h"


// This is an example of an exported variable
BASIC_API int nBasic=0;

// This is an example of an exported function.
BASIC_API int fnBasic(void)
{
	return 42;
}
