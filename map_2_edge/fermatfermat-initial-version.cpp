//Implements connected Fermat Sprial 
//Prints things on postscripts
//install ghost script to visualize the output https://www.ghostscript.com/

//We read the input and then scale it into a unit box 
//then scale everything to large int because this is how 
//clipperlib works 
//That way since the range of float number is bounded [0,1]
//the accuracy can be increased using large INT_SCALE (in def.h)


#include "fermat2.h"
#include <opencv2/highgui/highgui.hpp>
auto viz_fermat2 = new VIZ();

void Fermat::shape_cleanup(ClipperLib::Paths&p, const int tolerance){
	//clean up
	//some contours have points that are too close to one another
	//if so, we can remove one of them
	using namespace ClipperLib;
	for (int k = p[0].size() - 1; k >= 0; k--){
		for (int j = 0; j < p[0].size(); j++){
			if (j == k){ continue; }
			cInt d = dist(p[0][j].X, p[0][j].Y, cInt(0),
				p[0][k].X, p[0][k].Y, cInt(0));
			if (d <= tolerance){
				p[0].erase(p[0].begin() + k);
				break;
			}

		}
	}

}

void Fermat::contour_cleanup_too_close(std::vector<ClipperLib::Paths>&p, const int tolerance){
	//clean up
	//some contours have points that are too close to one another
	//if so, we can remove one of them
	using namespace ClipperLib;
	for (int i = 0; i < p.size(); i++){
		for (int j = 0; j < p[i].size(); j++){

			for (int k = p[i][j].size() - 1; k >= 1; k--){
				cInt d = dist(p[i][j][k].X, p[i][j][k].Y, cInt(0),
					p[i][j][k - 1].X, p[i][j][k - 1].Y, cInt(0));
				if (d <= tolerance){
					p[i][j].erase(p[i][j].begin() + k);
				}

			}
		}
	}
}

int Fermat::contour_path_check(const ClipperLib::Path&contour, std::vector<ClipperLib::Paths>paths){
	//we find the path on paths that this contour belongs to 
	//we do this by finding the closest point on each path in paths to contour 
	//we return the path with closest point 
	//we only check the outer path in paths (path[X].back())

	ClipperLib::cInt min_d = std::numeric_limits<ClipperLib::cInt>::max();;
	int min_d_path_id;
	for (int i = 0; i < paths.size(); i++){
		for (int j = 0; j < paths[i].back().size(); j++){
			for (int k = 0; k < contour.size(); k++){
				ClipperLib::cInt dd = dist(contour[k].X, contour[k].Y, ClipperLib::cInt(0),
					paths[i].back().at(j).X, paths[i].back().at(j).Y, ClipperLib::cInt(0));
				if (dd < min_d){
					min_d = dd;
					min_d_path_id = i;
				}
			}
		}
	}

	return min_d_path_id;
}

std::vector<ClipperLib::Paths> Fermat::reorder_path(const std::vector<ClipperLib::Paths>&p_sol){

	//create a new vector with Paths that seperate the parallel contour path

	std::vector<ClipperLib::Paths> p_contours;
	std::vector<ClipperLib::Paths> current(1);

	for (int i = 0; i < p_sol.size(); i++){
		std::cout << "p_sol.size()：" << p_sol.size() << std::endl;
		if (p_sol[i].size() == current.size()){

			std::cout << "p_sol[" << i << "].size()：" << p_sol[i].size() << " current.size()：" << current.size() << std::endl;

			if (p_sol[i].size() == 1){
				//one-to-one 
				current[0].push_back(p_sol[i][0]);
				std::cout << "current[0].push_back" << std::endl;
				std::cout << "*************" << std::endl;
			}
			else{
				for (int j = 0; j < p_sol[i].size(); j++){
					//check which (current) contour this path belongs to 
					int my_path_id = contour_path_check(p_sol[i][j], current);
					std::cout << "my_path_id：" << my_path_id << std::endl;
					current[my_path_id].push_back(p_sol[i][j]);

					//viz_concat_path(p_sol[i], true, true, false);					
					if (false){
						viz_fermat2->viz_concat_path(current[my_path_id], true, true, false);
					}
					if (false){
						viz_fermat2->path_plotter(p_sol[i][j], 0, 0, 0, 1, 1, 0, 1);
					}
				}
			}
		}
		else if (p_sol[i].size() >1 && current.size() == 1){
			//first add what you have to p_smart
			p_contours.push_back(current[0]);
			current.erase(current.begin(), current.end());
			std::cout << "p_contours.push_back(current[0]);" << std::endl;
			//add as many paths as you need 
			current.resize(p_sol[i].size());
			for (int j = 0; j < p_sol[i].size(); j++){
				current[j].push_back(p_sol[i][j]);
			}
		}
		else if (p_sol[i].size() < current.size() && current.size() > 1){

			//contruction, one of the current contours should be 
			//pushed into p_contours and new paths should be 
			//add to the remaining current path

			std::vector<bool> push(current.size(), true);

			for (int j = 0; j < p_sol[i].size(); j++){
				//check which (current) contour this path belongs to 
				int my_path_id = contour_path_check(p_sol[i][j], current);
				current[my_path_id].push_back(p_sol[i][j]);
				push[my_path_id] = false;

				//viz_concat_path(p_sol[i], true, true, false);					
				if (false){
					viz_fermat2->viz_concat_path(current[my_path_id], true, true, false);
				}
				if (false){
					viz_fermat2->path_plotter(p_sol[i][j], 0, 0, 0, 1, 1, 0, 1);
				}
			}

			for (int j = current.size() - 1; j >= 0; j--){
				if (push[j]){
					p_contours.push_back(current[j]);
					current.erase(current.begin() + j);
				}
			}

		}
		else{
			std::cout << "Error:: new case!!!" << std::endl;
			system("pause");
		}
	}

	for (int i = 0; i < current.size(); i++){
		p_contours.push_back(current[i]);
	}

	return p_contours;


}

void Fermat::contour_cleanup_too_short(std::vector<ClipperLib::Paths>&p_contours, const int threshold)
{
	//few paths have very few points. we remove such paths if the number 
	//of points is less than the threshold 
	for (int i = p_contours.size() - 1; i >= 0; i--){
		for (int j = p_contours[i].size() - 1; j >= 0; j--){
			if (p_contours[i][j].size() < threshold){
				p_contours[i].erase(p_contours[i].begin() + j);
			}
		}
		if (p_contours[i].size() < 2){
			p_contours.erase(p_contours.begin() + i);
		}
	}
}
void Fermat::init_fermat(std::vector<std::vector<cv::Point> > contours){
	srand(time(NULL));
	ClipperLib::Paths p;
	int clipper_offset_num;

	double  width_pixel_apart = 65.0;//pixel
	double  height_pixel_apart = 37.0;//pixel
	//std::vector<ClipperLib::Paths> p_sol(num_contours);

	ReadContour(p, clipper_offset_num,width_pixel_apart,height_pixel_apart, contours);

	ClipperLib::cInt width_pixel_apart_cInt = to_int<double, ClipperLib::cInt>(width_pixel_apart);
	ClipperLib::cInt height_pixel_apart_cInt = to_int<double, ClipperLib::cInt>(height_pixel_apart);

	// for (int pii = 0; pii < contours.size(); ++pii)
	// {
	//
	// 	ClipperLib::Path current_p;
	// 	current_p.clear();
	// 	p.clear();
	// 	for (int contour_j = 0; contour_j < contours[pii].size(); ++contour_j)
	// 	{
	// 		ClipperLib::IntPoint pt1;
	// 		//pt1.X = to_int<double, ClipperLib::cInt>(contours[contour_i][contour_j].x);//conver to ints 
	// 		//pt1.Y = to_int<double, ClipperLib::cInt>(contours[contour_i][contour_j].y);
	//
	// 		pt1.X = contours[pii][contour_j].x;//conver to ints 
	// 		pt1.Y = contours[pii][contour_j].y;
	// 		current_p.push_back(pt1);
	// 	}
	// 	p.push_back(current_p);
	// 	p_sol.push_back(p);
	// }

	//p_sol.clear();
	//for (int pii = 0; pii < p.size(); ++pii)
	//{

	//	ClipperLib::Paths temp_p;
	//	ClipperLib::Path current_p;
	//	temp_p.clear();
	//	current_p.clear();

	//	//for (int pij = 0; pij < p[pii].size; ++pij)
	//	//{

	//	//	ClipperLib::IntPoint pt1;
	//	//	//pt1.X = to_int<double, ClipperLib::cInt>(contours[contour_i][contour_j].x);//conver to ints 
	//	//	//pt1.Y = to_int<double, ClipperLib::cInt>(contours[contour_i][contour_j].y);

	//	//	pt1.X = p[pii][pij].X;//conver to ints 
	//	//	pt1.Y = p[pii][pij].Y;
	//	//	current_p.push_back(pt1);
	//	//}
	//	//p.push_back(current_p);
	//	std::cout << p[pii] << std::endl;
	//	p_sol.push_back(p);
	//}



	//path fill width 	
	//for airfoil_ex.poly, use 80 
	//for batman.poly, use 50
	ClipperLib::cInt max_x(-std::numeric_limits<ClipperLib::cInt>::max()),//当前系统最大值
		max_y(-std::numeric_limits<ClipperLib::cInt>::max()),
		min_x(std::numeric_limits<ClipperLib::cInt>::max()),
		min_y(std::numeric_limits<ClipperLib::cInt>::max());
	for (int p_0_size_i = 0; p_0_size_i < p[0].size(); ++p_0_size_i)
	{
		std::vector<ClipperLib::cInt> vec(2);
		vec[0] = (ClipperLib::cInt)p[0][p_0_size_i].X;
		vec[1] = (ClipperLib::cInt)p[0][p_0_size_i].Y;


		max_x = std::max(max_x, vec[0]);
		max_y = std::max(max_y, vec[1]);
		min_x = std::min(min_x, vec[0]);
		min_y = std::min(min_y, vec[1]);
	}



	//const ClipperLib::cInt clipper_offset_num_cInt = to_int<double, ClipperLib::cInt>(65.0);
	ClipperLib::cInt w;
	if (width_pixel_apart_cInt<height_pixel_apart_cInt)
	{
		w = width_pixel_apart_cInt;
	}else
	{
		w = height_pixel_apart_cInt;
	}

	//discretization spacing  (defined in the paper)	
	const int delta = 20000;
	const ClipperLib::cInt num_contours_clip = 20;

	viz_fermat2->viz_init("output5.ps", p);


	//compute the iso-contours 
	ClipperLib::ClipperOffset clip_offset;
	clip_offset.AddPaths(p, ClipperLib::jtRound, ClipperLib::etClosedPolygon);	
	std::vector<ClipperLib::Paths> p_sol(num_contours_clip);
	for (int i = 0; i < num_contours_clip; i++){
		clip_offset.Execute(p_sol[i], -(i+0.5)*w);//p_solution
	}
	for (int i = num_contours_clip - 1; i >= 0; i--){
		//squeeze up by removing path with zero segments 
		//(in case num_contours_clip is larger than neccessary)
		if (p_sol[i].size()==0){ 
			p_sol.pop_back(); 
		}
		else{ break; }
	}

	//clean up duplicated points on the contour path
	//contour_cleanup_too_close(p_sol, INT_TOL);


	//initialize the visualizer 
	//viz_init("output.ps", p);
	//add contours 
	//for (int i = 0; i < p_sol.size(); i++){
	//viz_concat_path(p_sol[i], false,true,false);
	//}

	//prepare vector to be passed to fermat function
	std::vector<ClipperLib::Paths> p_contours = reorder_path(p_sol);
	//for (int i = 0; i < p_contours.size(); i++){
	//	viz_concat_path(p_contours[i], true, true, false);
	//}

	//remove too short paths (mostly exists in the deep interior 
	//of a pocket)
	contour_cleanup_too_short(p_contours, 10);


	//viz_init("output.ps", p);
	//for (int i = 0; i < p_contours.size(); i++){
	//	viz_concat_path(p_contours[i], false, true, false);
	//}



	//reroute the contours into spirals
	std::vector<ClipperLib::Paths> fermat = Spirals::from_contour_to_fermat(p_contours,
		delta);





	return;
}

cv::Point Fermat::after_p(std::vector<cv::Point> p,cv::Point pt)
{
	//find the point before pt along the path p
	//pt should be on p
	int min_dist = std::numeric_limits<int>::max();
	int closest_i;
	for (int i = 0; i < p.size(); i++){
		int d = dist(int(pt.x), int(pt.y), int(0),
			int(p[i].x), int(p[i].y), int(0));
		if (d < min_dist){
			min_dist = d;
			closest_i = i;
		}
	}

	int bf = (closest_i == p.size()-1) ? 0 : closest_i + 1;//the one before the closest
	return p[bf];
}

int Fermat::find_index(std::vector<cv::Point> p, cv::Point pt)
{
	//find the index of pt along the path p
	//pt should be on p
	int min_dist = std::numeric_limits<int>::max();
	int closest_i;
	for (int i = 0; i < p.size(); i++){
		int d = dist(int(pt.x), int(pt.y), int(0),
			int(p[i].x), int(p[i].y), int(0));
		if (d < min_dist){
			min_dist = d;
			closest_i = i;
		}
	}
	return closest_i;
}

cv::Point Fermat::inside_p(std::vector<cv::Point> p,cv::Point pt)
{

	//find the point before pt along the path p
	//pt should be on p
	int min_dist = std::numeric_limits<int>::max();
	int closest_i;
	for (int i = 0; i < p.size(); i++){
		int d = dist(int(pt.x), int(pt.y), int(0),
			int(p[i].x), int(p[i].y), int(0));
		if (d < min_dist){
			min_dist = d;
			closest_i = i;
		}
	}
	return p[closest_i];
}
cv::Point Fermat::outside_p(std::vector<cv::Point> p,cv::Point pt)
{
	//find the point before pt along the path p
	//pt should be on p
	int min_dist = std::numeric_limits<int>::max();
	int closest_i;
	for (int i = 0; i < p.size(); i++){
		int d = dist(int(pt.x), int(pt.y), int(0),
			int(p[i].x), int(p[i].y), int(0));
		if (d < min_dist){
			min_dist = d;
			closest_i = i;
		}
	}
	return p[closest_i];
}

cv::Point Fermat::before_p(std::vector<cv::Point> p, cv::Point pt)
{
	//find the point before pt along the path p
	//pt should be on p
	int min_dist = std::numeric_limits<int>::max();
	int closest_i;
	for (int i = 0; i < p.size(); i++){
		int d = dist(int(pt.x), int(pt.y), int(0),
			int(p[i].x), int(p[i].y), int(0));
		if (d < min_dist){
			min_dist = d;
			closest_i = i;
		}
	}


	int bf = (closest_i == 0) ? p.size() - 1 : closest_i-1;//the one before the closest
	return p[bf];
}

void Fermat::get_fermat_from_contour(std::vector<std::vector<cv::Point> > contours)
{
	//2) convert spiral to fermat 

	int next_pt_id = 1;//id in the contour 


	int contourSize = contours.size();
	std::vector<std::vector<cv::Point> > contours_joint_point_first;
	std::vector<std::vector<cv::Point> > contours_joint_point_second;

	std::vector<cv::Point> result1;
	result1.clear();
	std::vector<cv::Point> result2;
	result2.clear();

	contours_joint_point_first.clear();
	cv::Point p_in;
	cv::Point after_p_in;
	cv::Point p_out;
	cv::Point p_in_second;
	cv::Point p_last_first;
	cv::Point p_last_second;
	int p_in_index;
	cv::Point contour_end_pt,inside_contour_end_pt,after_inside_contour_end_pt,before_p_in,before_before_p_in;
	
	for (int i = 0; i < contours.size()-1; ++i)
	{
		if (i == 0)
		{
			p_in = contours[0].front();
			p_out = contours[0][contours[0].size() - 2];
			contour_end_pt = before_p(contours[0],p_out);
		}
		
		p_in_index = find_index(contours[i],p_in);

		if (p_in_index == 0)
		{
			//in-end
			for (int j = p_in_index; j < contours[i].size(); ++j)
			{
				result1.push_back(contours[i][j]);
				if (contours[i][j] == contour_end_pt)
				{
					inside_contour_end_pt = inside_p(contours[i+1],contour_end_pt);
					//end-secondIn
					result1.push_back(inside_contour_end_pt);
					after_inside_contour_end_pt = after_p(contours[i+1],inside_contour_end_pt);
					result1.push_back(after_inside_contour_end_pt);

					p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
					before_p_in = before_p(contours[i+2],p_in);

					before_before_p_in = before_p(contours[i+2],before_p_in);

					contour_end_pt = before_p(contours[i+2],before_before_p_in);
					i++;
					break;
				}
			}			
		}
		else
		{
			//in-end
			for (int j = p_in_index; j < contours[i].size(); ++j)
			{
				result1.push_back(contours[i][j]);
				if (contours[i][j] == contour_end_pt)
				{
					inside_contour_end_pt = inside_p(contours[i+1],contour_end_pt);
					//end-secondIn
					result1.push_back(inside_contour_end_pt);
					after_inside_contour_end_pt = after_p(contours[i+1],inside_contour_end_pt);
					result1.push_back(after_inside_contour_end_pt);

					//secondIn-thirdIn
					if ((i+2) == contourSize-1)
					{
						p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
						contour_end_pt = before_p(contours[i+2],p_in);
					}	
					else
					{
						p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
						before_p_in = before_p(contours[i+2],p_in);

						before_before_p_in = before_p(contours[i+2],before_p_in);

						contour_end_pt = before_p(contours[i+2],before_before_p_in);
					}
					i++;
					break;
				}
			}
			//in-end
			for (int j = 0; j < p_in_index; ++j)
			{
				result1.push_back(contours[i][j]);
				if (contours[i][j] == contour_end_pt)
				{
					inside_contour_end_pt = inside_p(contours[i+1],contour_end_pt);
					//end-secondIn
					result1.push_back(inside_contour_end_pt);
					after_inside_contour_end_pt = after_p(contours[i+1],inside_contour_end_pt);
					result1.push_back(after_inside_contour_end_pt);

					//secondIn-thirdIn
					if ((i+2) == contourSize-1)
					{
						p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
						contour_end_pt = before_p(contours[i+2],p_in);
					}	
					else
					{
						p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
						before_p_in = before_p(contours[i+2],p_in);

						before_before_p_in = before_p(contours[i+2],before_p_in);

						contour_end_pt = before_p(contours[i+2],before_before_p_in);
					}
					i++;
					break;
				}
			}
		}
	}

	//get last contour point
	p_in_index = find_index(contours[contourSize-1],p_in);
	for (int j = p_in_index; j < contours[contourSize-1].size(); ++j)
	{
		result1.push_back(contours[contourSize-1][j]);
		if (contours[contourSize-1][j] == contour_end_pt)
		{
			break;
		}
	}
	for (int j = 0; j < p_in_index; ++j)
	{
		result1.push_back(contours[contourSize-1][j]);
		if (contours[contourSize-1][j] == contour_end_pt)
		{
			break;
		}
	}
	
	cv::Mat image(880, 1680, CV_8UC3, cv::Scalar(0, 0, 0));
	cv::Mat image2(880, 1680, CV_8UC3, cv::Scalar(0, 0, 0));

	for (int i = 0; i < result1.size()-1; ++i)
	{
		cv::Point start = result1[i];
		cv::Point end = result1[i+1];
		cv::line(image, start, end, cv::Scalar(0, 255, 255));
	}

	cv::imshow("image",image);
	cvWaitKey(0);

	//p_out and after p_out
	result2.push_back(p_out);
	result2.push_back(after_p(contours[0],p_out));

	//second
	for (int i = 1; i < contours.size()-1; ++i)
	{
		if (i == 1)
		{
			p_in = inside_p(contours[1],after_p(contours[0],p_out));
			contour_end_pt = before_p(contours[1],before_p(contours[1],before_p(contours[1],p_in)));
		}
		
		p_in_index = find_index(contours[i],p_in);

		if (p_in_index == 0)
		{
			//in-end
			for (int j = p_in_index; j < contours[i].size(); ++j)
			{
				result2.push_back(contours[i][j]);
				if (contours[i][j] == contour_end_pt)
				{
					inside_contour_end_pt = inside_p(contours[i+1],contour_end_pt);
					
					result2.push_back(inside_contour_end_pt);
					//end-secondIn
					after_inside_contour_end_pt = after_p(contours[i+1],inside_contour_end_pt);
					result2.push_back(after_inside_contour_end_pt);

					p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
					before_p_in = before_p(contours[i+2],p_in);

					before_before_p_in = before_p(contours[i+2],before_p_in);

					contour_end_pt = before_p(contours[i+2],before_before_p_in);
					i++;
					break;
				}
			}           
		}
		else
		{
			//in-end
			for (int j = p_in_index; j < contours[i].size(); ++j)
			{
				result2.push_back(contours[i][j]);
				if (contours[i][j] == contour_end_pt)
				{
					inside_contour_end_pt = inside_p(contours[i+1],contour_end_pt);
					result2.push_back(inside_contour_end_pt);

					//end-secondIn
					if ((i+2) >= contourSize)
					{					
						i++;
						break;
					}
					else
					{
						after_inside_contour_end_pt = after_p(contours[i+1],inside_contour_end_pt);
						result2.push_back(after_inside_contour_end_pt);

						p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
						before_p_in = before_p(contours[i+2],p_in);

						before_before_p_in = before_p(contours[i+2],before_p_in);

						contour_end_pt = before_p(contours[i+2],before_before_p_in);					
					}
					i++;
					break;
				}
			}
			//in-end
			for (int j = 0; j < p_in_index; ++j)
			{
				result2.push_back(contours[i][j]);
				if (contours[i][j] == contour_end_pt)
				{
					inside_contour_end_pt = inside_p(contours[i+1],contour_end_pt);
					result2.push_back(inside_contour_end_pt);

					//end-secondIn
					if ((i+2) >= contourSize)
					{					
						i++;
						break;
					}
					else
					{
						after_inside_contour_end_pt = after_p(contours[i+1],inside_contour_end_pt);
						result2.push_back(after_inside_contour_end_pt);

						p_in = inside_p(contours[i+2],after_inside_contour_end_pt);
						before_p_in = before_p(contours[i+2],p_in);

						before_before_p_in = before_p(contours[i+2],before_p_in);

						contour_end_pt = before_p(contours[i+2],before_before_p_in);					
					}
					i++;
					break;
				}
			}
		}
	}

	//for (int k = result2.size()-1; k >= 0; --k)
	//{
	//	cv::Point tmp = result2[k];
	//	result1.push_back(tmp);
	//}
	for (int i = 0; i < result2.size()-1; ++i)
	{
		cv::Point start = result2[i];
		cv::Point end = result2[i+1];
		cv::line(image, start, end, cv::Scalar(255, 255, 255));
	}

	cv::imshow("image2",image);
	cvWaitKey(0);

	// for (int i = 0; i < contours.size()-1; ++i)
	// {
	// 	std::vector<cv::Point> tmp;
	// 	tmp.clear();
 //
	// 	if (i==0)
	// 	{
	// 		p_in = contours[0].front();
	// 		p_in_second = contours[0][contours[0].size() - 2];
	// 	}
	// 	cv::Point p_in_before = before_p(contours[i],p_in);
	// 	cv::Point p_in_before_inside = inside_p(contours[i+1],p_in_before);
 //
	// 	tmp.push_back(p_in);
	// 	tmp.push_back(p_in_before);
 //
 //
	// 	if (i == contours.size()-2)
	// 	{
	// 		p_last_first = p_in_before_inside;
	// 	}
	// 	contours_joint_point_first.push_back(tmp);
	// 	p_in = p_in_before_inside;
	// }
 //
	// for (int i = 0; i < contours.size()-1; ++i)
	// {
	// 	if (i==0)
	// 	{
	// 		p_in = p_in_second;
	// 	}
	// 	cv::Point p_in_before = before_p(contours[i],p_in);
	// 	cv::Point p_in_before_inside = inside_p(contours[i+1],p_in_before);
	// 	std::vector<cv::Point> tmp;
	// 	tmp.clear();
	// 	tmp.push_back(p_in);
	// 	tmp.push_back(p_in_before);
 //
	// 	if (i == contours.size()-2)
	// 	{
	// 		p_last_second = p_in_before_inside;
	// 	}
	// 	contours_joint_point_second.push_back(tmp);
 //
	// 	p_in = p_in_before_inside;
	// }






	std::vector<int> reroute_pts;//index of the rereouting points in the inward path 


}