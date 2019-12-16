#pragma once
#include <stdio.h>
#include <cstdint>
#include <random>
#include <string>
#include <iostream>
#include <numeric>
#include <time.h>

#include "def.h"
#include "clipper.hpp"
#include "viz.h"
#include "convertor.h"
#include "spirals.h"
#include "geom.h"
#include "read_contour.h"
#include <opencv2/core/core.hpp>

class Fermat
{
public:
	static void shape_cleanup(ClipperLib::Paths&p, const int tolerance);
	static void contour_cleanup_too_close(std::vector<ClipperLib::Paths>&p, const int tolerance);
	static int contour_path_check(const ClipperLib::Path&contour, std::vector<ClipperLib::Paths>paths);
	static std::vector<ClipperLib::Paths> reorder_path(const std::vector<ClipperLib::Paths>&p_sol);
	static void contour_cleanup_too_short(std::vector<ClipperLib::Paths>&p_contours, const int threshold = 20);
	static void init_fermat(std::vector<std::vector<cv::Point> > contours);

	static cv::Point after_p(std::vector<cv::Point> p,cv::Point pt);
	static int find_index(std::vector<cv::Point> p,cv::Point pt);
	static cv::Point inside_p(std::vector<cv::Point> p,cv::Point pt);
	static cv::Point inside_p_closest(std::vector<cv::Point> p, cv::Point pt);
	static cv::Point outside_p(std::vector<cv::Point> p,cv::Point pt);
	static cv::Point before_p(std::vector<cv::Point> p,
	cv::Point pt);
	static void get_fermat_from_contour(std::vector<std::vector<cv::Point> > contours);
	static void get_drone_path(std::vector<cv::Point> result);
	static bool get_angle_bet_three_points(cv::Point P1,cv::Point P2,cv::Point P3);

};
