#include "viz.h"


double viz_xmin, viz_ymin, viz_xmax, viz_ymax, viz_scale;
std::string viz_file_name;

void VIZ::get_path_minmax(const ClipperLib::Path& p,
	double &xmin, double&ymin,
	double &xmax, double&ymax){

	xmin = std::numeric_limits<double>::max();
	ymin = std::numeric_limits<double>::max();
	xmax = -std::numeric_limits<double>::max();
	ymax = -std::numeric_limits<double>::max();

	for (int j = 0; j < p.size(); j++){
		double x(to_float<ClipperLib::cInt, double>(p[j].X)), y(to_float<ClipperLib::cInt, double>(p[j].Y));

		xmax = std::max(xmax, x);
		ymax = std::max(ymax, y);

		xmin = std::min(xmin, x);
		ymin = std::min(ymin, y);
	}

}

void VIZ::dot_plotter(const int x,
	const int y,
	const double r,
	const double g,
	const double b,
	bool is_big ){

	std::fstream file(viz_file_name, std::ios::app);

	double x1(to_float<ClipperLib::cInt, double>(x)), y1(to_float<ClipperLib::cInt, double>(y));

	file << r << "  " << g << "  " << b << "  "
		<< (x1 - viz_xmin)*viz_scale << " " << (y1 - viz_ymin)*viz_scale << " ";
	if (is_big){
		file << viz_scale / 200.0 << " color_dot" << std::endl;
	}
	else{
		file << viz_scale / 500.0 << " color_dot" << std::endl;
	}

	//file.close();
}

void VIZ::path_plotter(const ClipperLib::Path&p,
	const double r ,
	const double g ,
	const double b ,
	bool with_dots ,
	bool is_closed ,
	bool show_numbers ,
	bool is_thick ){

	std::fstream file(viz_file_name, std::ios::app);

	for (int j = 0; j < p.size() - 1; j++){
		double x1(to_float<ClipperLib::cInt, double>(p[j].X)), y1(to_float<ClipperLib::cInt, double>(p[j].Y));
		double x2(to_float<ClipperLib::cInt, double>(p[j + 1].X)), y2(to_float<ClipperLib::cInt, double>(p[j + 1].Y));

		file << SEG_WIDTH + SEG_WIDTH*is_thick << " "
			<< r << "  " << g << "  " << b << "  "
			<< (x1 - viz_xmin)*viz_scale << " " << (y1 - viz_ymin)*viz_scale << " "
			<< (x2 - viz_xmin)*viz_scale << " " << (y2 - viz_ymin)*viz_scale << " "
			<< "color_seg" << std::endl;

		x1 = x2;
		y1 = y2;
	}
	if (is_closed){

		double x1(to_float<ClipperLib::cInt, double>(p[0].X)), y1(to_float<ClipperLib::cInt, double>(p[0].Y));
		double x2(to_float<ClipperLib::cInt, double>(p[p.size() - 1].X)), y2(to_float<ClipperLib::cInt, double>(p[p.size() - 1].Y));

		file << SEG_WIDTH + SEG_WIDTH*is_thick << " "
			<< r << "  " << g << "  " << b << "  "
			<< (x1 - viz_xmin)*viz_scale << " " << (y1 - viz_ymin)*viz_scale << " "
			<< (x2 - viz_xmin)*viz_scale << " " << (y2 - viz_ymin)*viz_scale << " "
			<< "color_seg" << std::endl;
	}

	//file.close();

	if (with_dots){
		//plot all the dots at once 
		for (int j = 0; j < p.size(); j++){
			dot_plotter(p[j].X, p[j].Y, r, g, b);
		}
	}

	if (show_numbers){
		for (int j = 0; j < p.size(); j++){

			double x(to_float<ClipperLib::cInt, double>(p[j].X)),
				y(to_float<ClipperLib::cInt, double>(p[j].Y));

			file << "/Times-Roman findfont" << std::endl;
			file << "0.06 scalefont" << std::endl;
			file << "0 0 0 setrgbcolor" << std::endl;
			file << "setfont" << std::endl;
			file << (x - viz_xmin)*viz_scale << " " << (y - viz_ymin)*viz_scale << " moveto" << std::endl;
			//file << "(" << i<<"-"<<j<<") ";
			file << "(" << j << ") ";
			file << "show" << std::endl;
		}
	}
}

void VIZ::viz_init(const std::string file_name, const ClipperLib::Paths& p){

	std::fstream file(file_name, std::ios::out);
	viz_file_name = file_name;

#pragma region blocks 
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale     % one unit = one inch" << std::endl;

	file << "/quad_white      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " 1.0 1.0 1.0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_disk      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	//file << " fill" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/blue_disk      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0.33 1.0 setrgbcolor" << std::endl;
	//	file << " fill" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/purple_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.98 0 0.717 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.25 0.75 0.25 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/quad_purp      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0.49 0.15 0.8 setrgbcolor" << std::endl;
	//	file << "fill" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	//	file << "fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.006 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/color_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.006 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/black_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.006 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0 0 setrgbcolor" << std::endl;
	//	file <<47/255.0<<"	"<<79/255.0<<"	"<< 79/255.0<<" setrgbcolor" << std::endl;	
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << "0 0.392156862745098 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/purple_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.98 0 0.717 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/gray_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.83 0.83 0.83 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	//	file << " fill" << std::endl;
	file << " 0.003 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/r_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/faint_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 1.0 0.67 0.67 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/g_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/color_seg      % stack: x y r g b" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " setrgbcolor" << std::endl;
	file << " setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/blue_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 1 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "1 0 0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0.392156862745098 0 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;
#pragma endregion  

	viz_xmin = std::numeric_limits<double>::max();
	viz_ymin = std::numeric_limits<double>::max();
	viz_xmax = -std::numeric_limits<double>::max();
	viz_ymax = -std::numeric_limits<double>::max();

	for (int i = 0; i < p.size(); i++){
		double p_xmin, p_ymin, p_xmax, p_ymax;
		get_path_minmax(p[i], p_xmin, p_ymin, p_xmax, p_ymax);
		viz_xmin = std::min(viz_xmin, p_xmin);
		viz_ymin = std::min(viz_ymin, p_ymin);
		viz_xmax = std::max(viz_xmax, p_xmax);
		viz_ymax = std::max(viz_ymax, p_ymax);

	}

	double Lx(viz_xmax - viz_xmin), Ly(viz_ymax - viz_ymin), scale_x, scale_y, shift_x, shift_y;

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;
	//scale_x = 3.0 / Lx;
	//scale_y = 3.0 / Ly;

	if (scale_x < scale_y){
		viz_scale = scale_x;
		shift_x = 0.9;
		shift_y = (11.0 - (Lx*viz_scale)) / 2.0;
	}
	else{
		viz_scale = scale_y;
		//shift_x = (8.3-(Ly*scale))/2.0;
		shift_x = 0.9;
		shift_y = 1.35;
	}

	file << shift_x << " " << shift_y << " translate" << std::endl;
	file.close();

	for (int i = 0; i < p.size(); i++){
		path_plotter(p[i], 0.0, 0.0, 0.0, false, true, false, true);
	}

}

void VIZ::viz_concat_path(const ClipperLib::Paths& p,
	bool with_dots ,
	bool is_closed ,
	bool show_numbers,
	bool is_connected){
	//is_connected means the end point of one path is the 
	//starting point of the next one

	double r(double(rand()) / double(RAND_MAX)),
		g(double(rand()) / double(RAND_MAX)),
		b(double(rand()) / double(RAND_MAX));

	for (int i = 0; i < p.size(); i++){
		path_plotter(p[i], r, g, b, with_dots, is_closed, show_numbers);

		if (is_connected){
			if (i < int(p.size()) - 1){
				//connect the ends of inteior paths 
				std::fstream file(viz_file_name, std::ios::app);

				double x1(to_float<ClipperLib::cInt, double>(p[i].back().X)), y1(to_float<ClipperLib::cInt, double>(p[i].back().Y));
				double x2(to_float<ClipperLib::cInt, double>(p[i + 1].front().X)), y2(to_float<ClipperLib::cInt, double>(p[i + 1].front().Y));

				file << SEG_WIDTH << " "
					<< r << "  " << g << "  " << b << "  "
					<< (x1 - viz_xmin)*viz_scale << " " << (y1 - viz_ymin)*viz_scale << " "
					<< (x2 - viz_xmin)*viz_scale << " " << (y2 - viz_ymin)*viz_scale << " "
					<< "color_seg" << std::endl;

			}
		}
	}
}