#ifndef _READ_POLY_
#define _READ_POLY_

#include <fstream>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>

#include "def.h"
#include "convertor.h"

void ReadPoly(ClipperLib::Paths &p, std::string file_name){
	//read a single poly from file_name file and store it in clipper data 
	//structure Path
	//The input file format is .poly format 
	//ref: https://www.cs.cmu.edu/~quake/triangle.poly.html
	//Here we ignore  holes and only consider line segments 

	p.clear();

	std::ifstream input_file(file_name);
	if (!input_file.is_open()){
		PRINT_ERROR("Can not open the input file " + file_name);		
	}

	//**** Read geometry	
	int num_points, blah, is_boundary, id;	
	bool is_zero_index(true);

	input_file >> num_points >> blah >> blah >> is_boundary;//don't care about dim, attributes or boundary marker	
	std::vector<std::vector<double>> coords(num_points, std::vector<double>(2));	
	double max_x(-std::numeric_limits<double>::max()),//当前系统最大值
		   max_y(-std::numeric_limits<double>::max()),
		   min_x( std::numeric_limits<double>::max()),
		   min_y( std::numeric_limits<double>::max());

	for (int i = 0; i < num_points; i++){
		std::vector<double> vec(2);
		input_file >> id >> vec[0] >> vec[1];
		if (is_boundary){
			input_file >> blah;
		}	
		if (i == 0){
			is_zero_index = (id == i);
		}
		max_x = std::max(max_x, vec[0]);
		max_y = std::max(max_y, vec[1]);

		min_x = std::min(min_x, vec[0]);
		min_y = std::min(min_y, vec[1]);

		coords[i]=vec;
	}

	//shift everything to first quadrant (min_x,min_y)=(0,0)
	for (int i = 0; i < num_points; i++){
		coords[i][0] -= min_x;
		coords[i][1] -= min_y;
	}
	max_x -= min_x;
	max_y -= min_y;
	min_x = 0;
	min_y = 0;
	double scale = std::max(max_x, max_y);
	//scale inside the unit box 
	for (int i = 0; i < num_points; i++){
		coords[i][0] /= scale;
		coords[i][1] /= scale;

#if _DEBUG
		//check the int-to-float conversion 
		double i1 = to_float<ClipperLib::cInt, double>(to_int<double, ClipperLib::cInt>(coords[i][0]));
		double i2 = to_float<ClipperLib::cInt, double>(to_int<double, ClipperLib::cInt>(coords[i][1]));

		if ((abs(coords[i][0]) > 0.1 && abs(i1 - coords[i][0]) / coords[i][0]  > 0.1) ||
			(abs(coords[i][1]) > 0.1 && abs(i2 - coords[i][1]) / coords[i][1]  > 0.1)){
			PRINT_ERROR("Float-to-int conversion has large error");
		}

#endif
	}
	max_x /= scale;
	max_y /= scale;


	//**** Read topology
	int num_segments, start(std::numeric_limits<int>::max());
	input_file >> num_segments >> blah;	
	ClipperLib::Path current_p;
	int current_p_num_segments(0);
	for (int i =0; i < num_segments; i++){
		ClipperLib::cInt p1, p2;
		input_file >> blah >> p1 >> p2;
		if (is_boundary){
			input_file >> blah;
		}
		if (!is_zero_index){
			p1--;
			p2--;
		}
		if (start == std::numeric_limits<ClipperLib::cInt>::max()){
			//starting of a new path 
			start = p1;
			if (current_p.size() > 0){
				//store the previous path (if any) 
				p.push_back(current_p);
				ClipperLib::Path new_p;
				current_p.swap(new_p);
			}
		}
		else if (p2 == start){
			//end this path 
			start = std::numeric_limits<int>::max();
		}
		ClipperLib::IntPoint pt1;
		pt1.X = to_int<double, ClipperLib::cInt>(coords[p1][0]);//conver to ints 
		pt1.Y = to_int<double, ClipperLib::cInt>(coords[p1][1]);		
		current_p.push_back(pt1);
		
	}
	p.push_back(current_p);

}


#endif /*_READ_POLY_*/