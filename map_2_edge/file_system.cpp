
#include "file_system.h"

#include <iostream>
#include <fstream>



std::string FileSystem::extension(const std::string& path) {
	int len = int(path.length()) ;
	int point_pos  = -1 ;
	int base_start = -1 ;
	{for(int i=len-1; i>=0; i--) {
		if(point_pos == -1 && base_start == -1 && path[i] == '.') {
			point_pos = i ;
		}
		if(base_start == -1 && (path[i] == '/' || path[i] == '\\')) {
			base_start = i ;
		}
	}}
	if(point_pos == -1) {
		return "" ;
	}
	std::string result = path.substr(point_pos + 1, len - point_pos - 1) ;
	{for(unsigned int i=0; i<result.length(); i++) {
		result[i] = tolower(result[i]) ;
	}}
	return result ;
}

std::string FileSystem::base_name(const std::string& path) {
	int len = int(path.length()) ;
	int point_pos  = -1 ;
	int base_start = -1 ;
	for(int i=len-1; i>=0; i--) {
		if(point_pos == -1 && base_start == -1 && path[i] == '.') {
			point_pos = i ;
		}
		if(base_start == -1 && (path[i] == '/' || path[i] == '\\')) {
			base_start = i ;
		}
	}
	if(point_pos == -1) {
		point_pos = len ;
	}
	return path.substr(base_start + 1, point_pos - base_start - 1) ;
}

std::string FileSystem::dir_name(const std::string& path) {
	int len = int(path.length()) ;
	int base_start = -1 ;
	for(int i=len-1; i>=0; i--) {
		if(base_start == -1 && (path[i] == '/' ||path[i] == '\\')) {
			base_start = i ;
		}
	}
	if(base_start == -1) {
		return "." ;
	}
	std::string result = path.substr(0, base_start + 1) ;
	if(result[result.length() - 1] == '/') {
		result = result.substr(0, result.length() - 1) ;
	}
	return result ;
}

