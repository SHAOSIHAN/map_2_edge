#ifndef _GEOM_
#define _GEOM_

#include "def.h"
#include "convertor.h"

template <typename T>
void cross_prodcut(const T xv1,//input vector 1
	const T yv1,
	const T zv1,
	const T xv2,//input vector 2
	const T yv2,
	const T zv2,
	T&xx, T&yy, T&zz){//output vector 3

	//Find the cross product between vector 1 and vectro 2
	xx = yv1*zv2 - zv1*yv2;
	yy = zv1*xv2 - xv1*zv2;
	zz = xv1*yv2 - yv1*xv2;

}

template<typename T>
T dot_prodcut(const T xv1,
	const T yv1,
	const T zv1,
	const T xv2,
	const T yv2,
	const T zv2){

	//Dot product of two vectors 
	return xv1*xv2 + yv1*yv2 + zv1*zv2;
}

template <typename T>
T angleVecVec(const T xv1, const T yv1, const T zv1, const T xv2, const T yv2, const T zv2 ){

	T dot = dot_prodcut(xv1, yv1, zv1, xv2, yv2, zv2);
	if (dot == T(0)){ return 1.5707963268; }

	T angle = dot / sqrt((xv1*xv1 + yv1*yv1 + zv1*zv1) * (xv2*xv2 + yv2*yv2 + zv2*zv2));

	if (angle>1.0){ return 0.0; }
	if (angle<-1.0){ return PI; }

	return acos(angle);
}

template<typename T>
T dist(T x1, T y1, T z1, T x2, T y2, T z2)
{
	T dx, dy, dz;
	dx = x1 - x2;
	dy = y1 - y2;
	dz = z1 - z2;
	dx *= dx;
	dy *= dy;
	dz *= dz;

	return dx + dy + dz;
}
template <typename T>
inline bool point_line_segment_projection(const T ptx, const T pty,
	                                      const T x1, const T y1,
	                                      const T x2, const T y2,
	                                      T&xp, T&yp, T&distance){

	//return the point closest to (ptx, pty) on the line segment 
	//between (x1,y1) and (x2,y2)

	//scale everything to floats 
	//double ptx_d = to_float<int, double>(ptx);
	//double pty_d = to_float<int, double>(pty);

	//double x1_d = to_float<int, double>(x1);
	//double y1_d = to_float<int, double>(y1);

	//double x2_d = to_float<int, double>(x2);
	//double y2_d = to_float<int, double>(y2);

	//vector (x1,y1)->(x2,y2)
	T x12 = x2 - x1;
	T y12 = y2 - y1;

	//dist 
	T dd = dist<T>(x1, y1, 0, x2, y2, 0);

	if (abs(dd) < 1){
		PRINT_ERROR("line segment is actually a point!!!");
		return false;
		//(x1,y1) and (x2,y2) are the same point 
		xp = x1;
		yp = y1;
		distance = dist<T>(x1, y1, 0, ptx, pty, 0);
	}
	else{
		//vector (x1,y1)->(ptx,pty)
		T x1p = ptx - x1;
		T y1p = pty - y1;

		T t = dot_prodcut<T>(x1p, y1p, 0, x12, y12, 0);
		double t_d = double(t) / double(dd);

		if (t_d < -0.1 || t_d > 1.0 + 0.1){
			return false;
		}
		else{
			xp = x1 + t_d*x12;
			yp = y1 + t_d*y12;
			distance = dist<T>(xp, yp, 0, ptx, pty, 0);
			return true;
		}
	}
}

/*
inline bool point_line_segment_projection(const int ptx, const int pty,
	const int x1, const int y1,
	const int x2, const int y2,
	int&xp, int&yp, int&distance){

	//return the point closest to (ptx, pty) on the line segment
	//between (x1,y1) and (x2,y2)

	//scale everything to floats
	double ptx_d = to_float<int, double>(ptx);
	double pty_d = to_float<int, double>(pty);

	double x1_d = to_float<int, double>(x1);
	double y1_d = to_float<int, double>(y1);

	double x2_d = to_float<int, double>(x2);
	double y2_d = to_float<int, double>(y2);

	//vector (x1,y1)->(x2,y2)
	double x12 = x2_d - x1_d;
	double y12 = y2_d - y1_d;

	//dist
	double dd = dist(x1_d, y1_d, 0.0, x2_d, y2_d, 0.0);

	if (abs(dd) == TOL){
		PRINT_ERROR("line segment is actually a point!!!");
		return false;
		//(x1,y1) and (x2,y2) are the same point
		xp = x1;
		yp = y1;
		distance = to_int<double, int>(dist(x1_d, y1_d, 0.0, ptx_d, pty_d, 0.0));
	}
	else{
		//vector (x1,y1)->(ptx,pty)
		double x1p = ptx_d - x1_d;
		double y1p = pty_d - y1_d;

		double t = dot_prodcut(x1p, y1p, 0.0, x12, y12, 0.0);
		t /= dd;

		if (t < TOL || t > 1.0 - TOL){
			return false;
		}
		else{
			double xp_d = x1_d + t*x12;
			double yp_d = y1_d + t*y12;

			xp = to_int<double, int>(xp_d);
			yp = to_int<double, int>(yp_d);

			distance = dist(xp, yp, 0, ptx, pty, 0);
			return true;
		}
	}
}*/


#endif /*_GEOM_*/