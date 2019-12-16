#ifndef _VIZ_
#define _VIZ_

#include <fstream>
#include <string>
#include <limits>
#include <algorithm>
#include <random>
#include "clipper.hpp"
#include "convertor.h"
#include "geom.h"

#define SEG_WIDTH 0.008
class VIZ
{
public:

private:

public:
	static void get_path_minmax(const ClipperLib::Path& p,
		double &xmin, double&ymin,
		double &xmax, double&ymax);
	static void dot_plotter(const int x,
		const int y,
		const double r = 0.0,
		const double g = 0.0,
		const double b = 0.0,
		bool is_big = false);
	static void path_plotter(const ClipperLib::Path&p,
		const double r = 0.0,
		const double g = 0.0,
		const double b = 0.0,
		bool with_dots = true,
		bool is_closed = true,
		bool show_numbers = false,
		bool is_thick = false);
	static void viz_init(const std::string file_name, const ClipperLib::Paths& p);
	static void viz_concat_path(const ClipperLib::Paths& p,
		bool with_dots = true,
		bool is_closed = true,
		bool show_numbers = false,
		bool is_connected = false);
};

#endif /*_VIZ_*/