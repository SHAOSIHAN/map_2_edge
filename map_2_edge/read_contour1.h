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

inline void ReadContour(ClipperLib::Paths &p,std::vector<ClipperLib::Paths> &p_sol,int &clipper_offset_num, std::vector<std::vector<cv::Point> > contours)
{
	ClipperLib::Paths p_temp;
	ClipperLib::Path current_p;
	p_temp.clear();
	p.clear();
	if (!(contours.size()>0)){
		PRINT_ERROR("There is no contour!");
	}


	double max_x(-std::numeric_limits<double>::max()),//当前系统最大值
		max_y(-std::numeric_limits<double>::max()),
		min_x(std::numeric_limits<double>::max()),
		min_y(std::numeric_limits<double>::max()),		
		min_x_ori(std::numeric_limits<double>::max()),
		min_y_ori(std::numeric_limits<double>::max());
	double scale;

	for (int contour_out_i = 0; contour_out_i < contours.size(); ++contour_out_i)
	{
		p_temp.clear();
		current_p.clear();

		if (contour_out_i == 0)
		{
			int num_points = contours[contour_out_i].size();
			std::vector<std::vector<double>> coords(num_points, std::vector<double>(2));

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
			min_x_ori = min_x;
			min_y_ori = min_y;
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
			scale = std::max(max_x, max_y);
			//scale inside the unit box 
			
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
		}
		else
		{
			int num_points = contours[contour_out_i].size();
			std::vector<std::vector<double>> coords(num_points, std::vector<double>(2));

			for (int contour_in_j = 0; contour_in_j < num_points; ++contour_in_j)
			{
				std::vector<double> vec(2);
				vec[0] = (double)contours[contour_out_i][contour_in_j].x;
				vec[1] = (double)contours[contour_out_i][contour_in_j].y;

				coords[contour_in_j] = vec;
			}


			//shift everything to first quadrant (min_x,min_y)=(0,0)
			for (int i = 0; i < num_points; i++){
				coords[i][0] -= min_x_ori;
				coords[i][1] -= min_y_ori;
			}

			//scale inside the unit box 
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
		}

		p_temp.push_back(current_p);
		p.push_back(current_p);
		//break;
		p_sol.push_back(p_temp);
	}
}

#endif