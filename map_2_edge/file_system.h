#ifndef _FILE_SYSTEM_H_
#define _FILE_SYSTEM_H_

//#include "basic_common.h"
#include <string>
#include <vector>


class FileSystem 
{
public:
	static std::string extension(const std::string& path) ;
	static std::string base_name(const std::string& path) ;
	static std::string dir_name(const std::string& path) ;
} ;


#endif
