#ifndef _READ_CONTOUR_
#define _READ_CONTOUR_

#include <fstream>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include "clipper.hpp"
#include "def.h"
#include <opencv2/core/core.hpp>

inline void ReadContour(ClipperLib::Paths &p, int &clipper_offset_num,double  &width_pixel_apart,double  &height_pixel_apart, std::vector<std::vector<cv::Point> > contours)
{
	p.clear();

	if (!(contours.size() > 0)){
		PRINT_ERROR("There is no contour!");
	}

	for (int contour_out_i = 0; contour_out_i < contours.size(); ++contour_out_i)
	{
		int num_points = contours[contour_out_i].size();
		std::vector<std::vector<double>> coords(num_points, std::vector<double>(2));
		double max_x(-std::numeric_limits<double>::max()),//ǰϵͳֵ
			max_y(-std::numeric_limits<double>::max()),
			min_x(std::numeric_limits<double>::max()),
			min_y(std::numeric_limits<double>::max());
		for (int contour_in_j = 0; contour_in_j < num_points; ++contour_in_j)
		{
			std::vector<double> vec(2);
			vec[0] = (double)contours[contour_out_i][contour_in_j].x;
			vec[1] = (double)contours[contour_out_i][contour_in_j].y;


			max_x = std::max(max_x, vec[0]);
			max_y = std::max(max_y, vec[1]);
			min_x = std::min(min_x, vec[0]);
			min_y = std::min(min_y, vec[1]);

			coords[contour_in_j] = vec;
		}

		//use max to compute fill width
		clipper_offset_num = (max_x - min_x) / 65;

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
		ClipperLib::Path current_p;
		width_pixel_apart = width_pixel_apart / scale;
		height_pixel_apart = height_pixel_apart / scale;
		for (int i = 0; i < num_points; i++){
			coords[i][0] /= scale;
			coords[i][1] /= scale;

#if _DEBUG
			//check the int-to-float conversion 
			double i1 = to_float<ClipperLib::cInt, double>(to_int<double, ClipperLib::cInt>(coords[i][0]));
			double i2 = to_float<ClipperLib::cInt, double>(to_int<double, ClipperLib::cInt>(coords[i][1]));

			if ((abs(coords[i][0]) > 0.1 && abs(i1 - coords[i][0]) / coords[i][0]  > 0.1) ||
				(abs(coords[i][1]) > 0.1 && abs(i2 - coords[i][1]) / coords[i][1] > 0.1)){
				PRINT_ERROR("Float-to-int conversion has large error");
			}

#endif

			ClipperLib::IntPoint pt1;
			pt1.X = to_int<double, ClipperLib::cInt>(coords[i][0]);//conver to ints 
			pt1.Y = to_int<double, ClipperLib::cInt>(coords[i][1]);
			current_p.push_back(pt1);

		}
		max_x /= scale;
		max_y /= scale;


		p.push_back(current_p);
		break;
	}
}

#endif