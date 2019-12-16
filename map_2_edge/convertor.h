#ifndef _CONVERTOR_
#define _CONVERTOR_
#include <string>

#include "def.h"

template <typename T_in, typename T_out>
//example: T_in = double , T_out = int
T_out to_int(T_in val){
	//scale double values to int 
	//intT could be int or long int or long long int 

	if (val < 0){		
		T_out  new_val = static_cast<T_out>(val*INT_SCALE - 0.5);

		/*if (new_val < std::numeric_limits<T_out>::min() ||
			new_val > std::numeric_limits<T_out>::max()){
			PRINT_ERROR(" Over/underflow. The value you want to convert it too high");
		}*/
		return new_val;
	}
	else{
		T_out  new_val = static_cast<T_out>(val*INT_SCALE + 0.5);

		/*if (new_val < std::numeric_limits<T_out>::min() ||
			new_val > std::numeric_limits<T_out>::max()){
			PRINT_ERROR(" Over/underflow. The value you want to convert it too high");
		}*/
		return new_val;
	}
}

template <typename T_in, typename T_out>
//example: T_in = int , T_out = double
T_out to_float(T_in val){
	//scale the int back to double 

	if (val < 0){
		T_out new_val = static_cast<T_out> ((val + 0.5) *FLOAT_SCALE);
		return new_val;
	}
	else{
		T_out new_val = static_cast<T_out>((val - 0.5) *FLOAT_SCALE);
		return new_val;
	}
}

#endif /*_CONVERTOR_*/