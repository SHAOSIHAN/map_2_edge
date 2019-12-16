#ifndef _SPIRAL_REROUTE_
#define _SPIRAL_REROUTE_

#include "clipper.hpp"
#include "convertor.h"
#include "geom.h"
#include "viz.h"
#include <iostream>

#define GAP 0.8

class Spirals
{
public:
	static bool point_contour_projection(ClipperLib::IntPoint pt,
		ClipperLib::Path p,
		int&seg_st, int&seg_end,
		ClipperLib::IntPoint&projection,
		bool compute_closest_point);
	static ClipperLib::IntPoint before_pt(const ClipperLib::Path&p,
		const ClipperLib::IntPoint pt,
		const int delta);
	static std::vector<ClipperLib::Paths> get_spiral(const std::vector<ClipperLib::Paths>&p, const int delta);

	static ClipperLib::Paths get_fermat(const ClipperLib::Paths&spiral,
		const int delta);
	static std::vector<ClipperLib::Paths> from_contour_to_fermat(
		const std::vector<ClipperLib::Paths>&p, const int delta = 20000);
protected:
private:
};


#endif /*_SPIRAL_REROUTE_*/