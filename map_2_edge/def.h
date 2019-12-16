//*** Parameters def
#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

#ifndef TwoPI
#define TwoPI 6.2831853071795862
#endif 

#ifndef FLOAT_TOL
#define FLOAT_TOL 10E-8
#endif 

#ifndef INT_TOL
#define INT_TOL 10
#endif 

#ifndef DIM
#define DIM 2
#endif

#ifndef RadToDeg
#define RadToDeg 57.295779513078550 //180.0/PI (muliply this by the radian angle to convert to degree)
#endif 

#ifndef INT_SCALE
//used for conversion from double to int
#define INT_SCALE (2<<12)
#endif 

#ifndef FLOAT_SCALE
#define FLOAT_SCALE (1.0/INT_SCALE)
#endif
//*********************

//*** PRINT_ERROR
#ifndef PRINT_ERROR
#include <stdio.h>
#include <string>
inline void Err(std::string err_line, const char *file, int line) {
	//Display the err_line 	
	printf("Error::%s \n Error generated in %s at line %d\n", err_line.c_str(), file, line);

#ifdef _WIN32
	system("pause");
#else
	exit(EXIT_FAILURE);
#endif

}
#define PRINT_ERROR( err_line ) (Err( err_line, __FILE__, __LINE__ ))
#endif 
//*********************


