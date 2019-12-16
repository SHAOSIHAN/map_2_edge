#include "spirals.h"

bool Spirals::point_contour_projection(ClipperLib::IntPoint pt,
	ClipperLib::Path p,
	int&seg_st, int&seg_end,
	ClipperLib::IntPoint&projection,
	bool compute_closest_point){
	//project pt into the countour path p
	//using brute force i.e., check against all line segments 
	//if compute_closest_point is set to true, then the method will
	//return the mid point of closest segment in case orthogonal projection failed 
	//TODO accelerate the computation here 


	ClipperLib::cInt min_dist_projection = std::numeric_limits<int>::max();
	ClipperLib::cInt min_dist_point = std::numeric_limits<int>::max();
	int min_dist_point_id;

	int seg_end_id = p.size() - 1;
	for (int seg_start_id = 0; seg_start_id < p.size(); seg_start_id++){

		ClipperLib::IntPoint segment_start;
		segment_start.X = p[seg_start_id].X;
		segment_start.Y = p[seg_start_id].Y;

		ClipperLib::IntPoint segment_end;
		segment_end.X = p[seg_end_id].X;
		segment_end.Y = p[seg_end_id].Y;

		ClipperLib::cInt dd_pt = dist(pt.X, pt.Y, ClipperLib::cInt(0),
			segment_start.X, segment_start.Y, ClipperLib::cInt(0));
		if (dd_pt < min_dist_point){ min_dist_point = dd_pt; min_dist_point_id = seg_start_id; }
		dd_pt = dist(pt.X, pt.Y, ClipperLib::cInt(0),
			segment_end.X, segment_end.Y, ClipperLib::cInt(0));
		if (dd_pt < min_dist_point){ min_dist_point = dd_pt; min_dist_point_id = seg_end_id; }

		ClipperLib::cInt xp, yp, dd;
		if (point_line_segment_projection(pt.X, pt.Y,
			segment_start.X, segment_start.Y,
			segment_end.X, segment_end.Y,
			xp, yp, dd)){
			if (dd < min_dist_projection){
				min_dist_projection = dd;
				projection.X = xp;
				projection.Y = yp;
				seg_st = seg_start_id;
				seg_end = seg_end_id;
			}
		}

		/*		std::cout << "point_contour_projection" <<" pt£º"
				<<pt.X<<" "<<pt.Y
				<<" xp: "<<xp
				<<" yp: "<<yp << std::endl;	*/
		seg_end_id = seg_start_id;
	}

	if (compute_closest_point){
		//if orth projection did not work 
		//then just connect to the closest point Z
		//Z is mid point between min_dist_point_id and second closest 
		//point to pt that is connect to min_dist_point_id
		int min_dist_point_id_next, min_dist_point_id_before;
		if (min_dist_point_id == 0){
			//it is only one segment 0 <-> 1 
			min_dist_point_id_next = min_dist_point_id_before = 1;
		}
		else if (min_dist_point_id == p.size() - 1){
			//it is only one segment [p.size() - 1] <-> [p.size() - 1]
			min_dist_point_id_next = min_dist_point_id_before = p.size() - 2;
		}
		else{
			min_dist_point_id_next = (min_dist_point_id == p.size() - 1) ? 0 : min_dist_point_id + 1;
			min_dist_point_id_before = (min_dist_point_id == 0) ? p.size() - 1 : min_dist_point_id - 1;
		}

		ClipperLib::cInt dd_next = dist(pt.X, pt.Y, ClipperLib::cInt(0),
			p[min_dist_point_id_next].X,
			p[min_dist_point_id_next].Y, ClipperLib::cInt(0));

		ClipperLib::cInt dd_before = dist(pt.X, pt.Y, ClipperLib::cInt(0),
			p[min_dist_point_id_before].X,
			p[min_dist_point_id_before].Y, ClipperLib::cInt(0));

		int second_closest = (dd_next < dd_before) ? min_dist_point_id_next : min_dist_point_id_before;
		ClipperLib::IntPoint close_pt;
		close_pt.X = (double(p[min_dist_point_id].X) + double(p[second_closest].X)) / 2.0;
		close_pt.Y = (double(p[min_dist_point_id].Y) + double(p[second_closest].Y)) / 2.0;
		//seg_st = min_dist_point_id;
		//seg_end = second_closest;		
		ClipperLib::cInt min_dist_closest = dist(pt.X, pt.Y, ClipperLib::cInt(0),
			close_pt.X, close_pt.Y, ClipperLib::cInt(0));

		if (min_dist_projection == std::numeric_limits<int>::max()){
			//if we have not any orth projection, then we return right away the 
			projection = close_pt;
			seg_st = min_dist_point_id;
			seg_end = second_closest;
		}
		else{
			//otherwise we check on the distance to recover from mistakes sometimes happens
			//we return the closest point only if the projection point is too far 
			//i.e., two times far 			
			if (min_dist_projection > 2 * min_dist_closest){
				projection = close_pt;
				seg_st = min_dist_point_id;
				seg_end = second_closest;
				return true;
			}
			else{
				//either projection is closer than the computed closest point 
				return true;
			}


		}
	}
	else {
		if (min_dist_projection == std::numeric_limits<int>::max()){
			return false;
		}
		else{
			return true;
		}
	}




}
ClipperLib::IntPoint Spirals::before_pt(const ClipperLib::Path&p,
	const ClipperLib::IntPoint pt,
	const int delta){
	//find the point before pt along the path p
	//pt should be on p
	ClipperLib::cInt min_dist = std::numeric_limits<int>::max();
	int closest_i;
	for (int i = 0; i < p.size(); i++){
		ClipperLib::cInt d = dist(ClipperLib::cInt(pt.X), ClipperLib::cInt(pt.Y), ClipperLib::cInt(0),
			ClipperLib::cInt(p[i].X), ClipperLib::cInt(p[i].Y), ClipperLib::cInt(0));
		if (d < min_dist){
			min_dist = d;
			closest_i = i;
		}
	}

	//move a distance delta from the closest point 
	int bf1 = closest_i;
	int bf = (closest_i == 0) ? p.size() - 1 : closest_i - 1;//the one before the closest
	ClipperLib::cInt distance_walked(0);
	while (true){
		ClipperLib::cInt d = dist(ClipperLib::cInt(p[bf].X), ClipperLib::cInt(p[bf].Y), ClipperLib::cInt(0),
			ClipperLib::cInt(p[bf1].X), ClipperLib::cInt(p[bf1].Y), ClipperLib::cInt(0));

		if (distance_walked + d >= delta){
			//interpolate and return the point 
			ClipperLib::cInt remaning_d = delta - distance_walked;

			ClipperLib::IntPoint return_pt;
			return_pt.X = p[bf1].X + double(double(remaning_d) / double(d))*(p[bf].X - p[bf1].X);
			return_pt.Y = p[bf1].Y + double(double(remaning_d) / double(d))*(p[bf].Y - p[bf1].Y);

			return return_pt;

		}
		else{
			distance_walked += d;
		}

		bf1 = bf;
		bf = (bf == 0) ? p.size() - 1 : bf - 1;
	}




}
std::vector<ClipperLib::Paths> Spirals::get_spiral(const std::vector<ClipperLib::Paths>&p, const int delta){

	//create spirals for each spirallable region in p	
	std::vector<ClipperLib::Paths> spiral(p.size());
	//for (int i = 0; i < p.size(); i++){
	int i = 0;
	while (i < p.size()){
		bool RE_DO = false;
		for (int j = 0; j < p[i].size(); j++){

			//copy all points to the spiral 
			//then project the last point to the next contour  


			//choose the starting point here of this contour 			
			ClipperLib::cInt k_start;

			if (j == 0){
				//could be picked randomly 
				//float rnd = (float(rand()) / float(RAND_MAX));
				//k_start = float(p[i][j].size() - 1) * rnd;

				//we need to pick a starting point such that it can be 
				//projected to all subsequent contours 
				//this is not perfect solution since for subsequent contours
				//we shift a little bit then the project points might
				//have no intersection and the spiral won't be completed 
				while (true){
					float rnd = (float(rand()) / float(RAND_MAX));
					std::cout << "rnd£º" << rnd << std::endl;
					k_start = float(p[i][j].size() - 1) * rnd;
					std::cout << "k_start£º" << k_start << std::endl;
					//batman
					//if (i == 0){ k_start = 1; break; }
					//else if (i == 1){ k_start = 33; break; }
					//else if (i == 2){ k_start = 39; break; }
					//else if (i == 3){ k_start = 77; break; }

					//if (i == 0){ k_start = 10; break; }
					//else if (i == 1){ k_start = 36; break; }
					//else if (i == 2){ k_start = 28; break; }



					ClipperLib::cInt k_before = (k_start == 0) ? p[i][j].size() - 1 : k_start - 1;
					ClipperLib::cInt k_after = (k_start == p[i][j].size() - 1) ? 0 : k_start + 1;

					//get the angle at this point (k_start)
					//if it is a sharp angle, then it is not suitable 
					//we need starting point at smooth part of the contour 
					double ang = RadToDeg * angleVecVec(double(p[i][j][k_before].X) - double(p[i][j][k_start].X),
						double(p[i][j][k_before].Y) - double(p[i][j][k_start].Y),
						0.0,
						double(p[i][j][k_after].X) - double(p[i][j][k_start].X),
						double(p[i][j][k_after].Y) - double(p[i][j][k_start].Y),
						0.0);
					if (ang < 100.0){
						continue;
					}


					bool has_projection_on_all_contours = true;
					for (int l = 1; l < p[i].size(); l++){
						std::cout << "l: " << l <<" of p["<<i<<"].size(): "<<p[i].size()<< std::endl;

						int seg_st, seg_end;
						ClipperLib::IntPoint projection;
						if (!point_contour_projection(p[i][j][k_start], p[i][l], seg_st, seg_end, projection, false)){
							std::cout << "projection" <<projection.X<<" "<<projection.Y << std::endl;
							has_projection_on_all_contours = false;
							break;
						}
					}
					if (has_projection_on_all_contours){ break; }

				}

			}
			else {
				//next contours starts from the projection of previous contour 
				//on this contour. we try to find the closest point on this
				//contour to the projection point such that it is on the 
				//right "side" 

				//from the previous contour, the last point and its projection
				//form a plane such that the first point on this previous contour 
				//is on the +ve side of this plane 
				//we want the closest point to the projected point on this new 
				//contour such that it falls in the +ve side of this plane 

				ClipperLib::IntPoint last;
				last = spiral[i][j - 1].at(spiral[i][j - 1].size() - 2);

				ClipperLib::IntPoint proj;
				proj = spiral[i][j - 1].back();

				ClipperLib::IntPoint first;
				first = spiral[i][j - 1].front();

				//https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
				ClipperLib::cInt first_sgn = (proj.X - last.X) * (first.Y - last.Y) -
					(proj.Y - last.Y) * (first.X - last.X);

				double min_dist = std::numeric_limits<int>::max();
				int min_k_id;
				for (int l = 0; l < p[i][j].size(); l++){
					ClipperLib::IntPoint pnt = p[i][j][l];//test point 

					ClipperLib::cInt p_sgn = (proj.X - last.X) * (pnt.Y - last.Y) -
						(proj.Y - last.Y) * (pnt.X - last.X);
					if (p_sgn *first_sgn > 0.0){
						//on the same side 
						ClipperLib::cInt dd = dist(ClipperLib::cInt(pnt.X), ClipperLib::cInt(pnt.Y), ClipperLib::cInt(0),
							ClipperLib::cInt(proj.X), ClipperLib::cInt(proj.Y), ClipperLib::cInt(0));
						if (dd < min_dist){
							min_dist = dd;
							min_k_id = l;
						}
					}
				}
				k_start = min_k_id;
			}


			ClipperLib::Path contour;
			ClipperLib::cInt k = k_start;
			ClipperLib::cInt k_end = (k == 0) ? p[i][j].size() - 1 : k - 1;
			ClipperLib::cInt k_plus;

			while (true){
				k_plus = (k == (p[i][j].size()) - 1) ? 0 : k + 1;

				ClipperLib::IntPoint pt;
				pt.X = (p[i][j][k].X + p[i][j][k_plus].X) / 2;
				pt.Y = (p[i][j][k].Y + p[i][j][k_plus].Y) / 2;

				if (j == p[i].size() - 1 && contour.size() > 20){
					//for the last contour, make sure there is enough clearance between 
					//the last point and the project point on this last contour
					//we wait until we add at least 20 point to activate this check 
					ClipperLib::cInt dddd = dist(pt.X, pt.X, ClipperLib::cInt(0),
						spiral[i].back().back().X, spiral[i].back().back().Y, ClipperLib::cInt(0));
					if (dddd < delta){
						break;
					}
				}
				contour.push_back(pt);

				k = k_plus;
				if ((j != 0 && k_plus == k_end) || (j == 0 && k_plus == k_start)){ break; }
			}

			if (false){
				VIZ::path_plotter(contour, 0, 0, 0, 1, 0, 0);
			}



			if (j < p[i].size() - 1){
				if (j > 0 && j < int(p[i].size()) - int(2)){
					//add more point to close large gaps 
					//ClipperLib::cInt dd = dist(contour.back().X, contour.back().Y, ClipperLib::cInt(0),
					//	                       spiral.back().back().X, spiral.back().back().Y, ClipperLib::cInt(0));

					ClipperLib::IntPoint pt;
					pt.X = contour.back().X + GAP*double(spiral[i].back().back().X - contour.back().X);
					pt.Y = contour.back().Y + GAP*double(spiral[i].back().back().Y - contour.back().Y);
					contour.push_back(pt);
				}

				//there is nothing to project at for the last contour 
				int seg_st, seg_end;
				ClipperLib::IntPoint projection;
				if (!point_contour_projection(contour.back(), p[i][j + 1], seg_st, seg_end, projection, false)){
					//because of choosing wrong starting point k_start, we can find good projection here
					//for that, we re-do everything for this spiral and wish for better luck 
					RE_DO = true;
					std::cout << " Re-doing spiral[" << i << "] " << std::endl;
					break;
					//PRINT_ERROR(" get_spiral():: Can not find projection. Spiral [" + std::to_string(i) + "] won'b be completed!!");
				}
				else {
					contour.push_back(projection);
				}

			}

			spiral[i].push_back(contour);

			if (false){
				VIZ::viz_concat_path(spiral[i], true, false, false, true);
			}
		}

		if (!RE_DO){ i++; }
		else {
			spiral[i].erase(spiral[i].begin(), spiral[i].end());

		}
	}

	return spiral;

}

ClipperLib::Paths Spirals::get_fermat(const ClipperLib::Paths&spiral,
	const int delta){

	//2) convert spiral to fermat 


	//the fermat will consist of two paths 
	//inward path and outward path
	ClipperLib::Path inward, outward;

	int current_path_id = 0;
	ClipperLib::Path current_path = spiral[current_path_id];

	int next_pt_id = 1;//id in the contour 

	ClipperLib::IntPoint p_in = spiral[0].front();
	ClipperLib::IntPoint p_out = spiral[0][spiral[0].size() - 2];

	ClipperLib::IntPoint next_pt;//next point to add to the fermat path
	next_pt = p_in;


	ClipperLib::IntPoint contour_end_pt;
	//first-contour end pt is the before(p_out)
	contour_end_pt = before_pt(spiral[current_path_id], p_out, delta);


	std::vector<int> reroute_pts;//index of the rereouting points in the inward path 
	//used to construct the outward link 

	inward.push_back(next_pt);
	int points_take_from_that_path = 0;
	while (true){
		next_pt = current_path[next_pt_id];
		points_take_from_that_path++;
		//check the distance to the before(contour_end_pt)
		ClipperLib::cInt d = dist(ClipperLib::cInt(contour_end_pt.X),
			ClipperLib::cInt(contour_end_pt.Y),
			ClipperLib::cInt(0),
			ClipperLib::cInt(next_pt.X),
			ClipperLib::cInt(next_pt.Y),
			ClipperLib::cInt(0));
		if (d < delta &&  int(current_path.size()) - points_take_from_that_path < 3){
			//here we add next_candidate and contour_end_pt
			inward.push_back(next_pt);
			inward.push_back(contour_end_pt);

			if (current_path_id > 0){
				reroute_pts.push_back(inward.size() - 1);
			}

			//increment to the next path
			if (current_path_id != spiral.size() - 1){
				current_path_id++;
				ClipperLib::Path next_path = spiral[current_path_id];
				next_path.swap(current_path);
				points_take_from_that_path = 0;
			}

			//then reroute and prepare the next contour 
			//by projecting to the next path 			
			int seg_st, seg_end;

			ClipperLib::IntPoint end_proj;
			if (!point_contour_projection(contour_end_pt, current_path, seg_st, seg_end, end_proj, true)){
				PRINT_ERROR("get_fermat():: Can not find the projection in the inward link(1)");
			}
			else{
				inward.push_back(end_proj);
			}



			int dd_end = dist(end_proj.X, end_proj.Y, ClipperLib::cInt(0),
				current_path[seg_end].X, current_path[seg_end].Y, ClipperLib::cInt(0));
			int dd_st = dist(end_proj.X, end_proj.Y, ClipperLib::cInt(0),
				current_path[seg_st].X, current_path[seg_st].Y, ClipperLib::cInt(0));

			if (dd_end < INT_TOL){
				next_pt_id = seg_end;
			}
			else if (dd_st < INT_TOL){
				next_pt_id = seg_st;
			}
			else{

				ClipperLib::cInt sgn = (contour_end_pt.X - next_pt.X) * (end_proj.Y - next_pt.Y) -
					(contour_end_pt.Y - next_pt.Y) * (end_proj.X - next_pt.X);
				ClipperLib::cInt sgn_st = (contour_end_pt.X - current_path[seg_st].X) * (end_proj.Y - current_path[seg_st].Y) -
					(contour_end_pt.Y - current_path[seg_st].Y) * (end_proj.X - current_path[seg_st].X);
				ClipperLib::cInt sgn_end = (contour_end_pt.X - current_path[seg_end].X) * (end_proj.Y - current_path[seg_end].Y) -
					(contour_end_pt.Y - current_path[seg_end].Y) * (end_proj.X - current_path[seg_end].X);
				//update next_pt_id with the one on the opposite of next_pt (so we can keep walking)
				if (sgn*sgn_st < 0 && sgn*sgn_end > 0){
					next_pt_id = seg_st;
				}
				else if (sgn*sgn_end < 0 && sgn*sgn_st > 0){
					next_pt_id = seg_end;
				}
				else{
					PRINT_ERROR("projected point is not on the right segment");
				}
			}
			//update the contour_end_pt 
			if (current_path_id < spiral.size() - 2){
				contour_end_pt = before_pt(spiral[current_path_id], end_proj, delta);

				if (!point_contour_projection(contour_end_pt, spiral[current_path_id + 1], seg_st, seg_end, contour_end_pt, true)){
					PRINT_ERROR("get_fermat():: Can not find the projection in the inward link(2)");
				}
			}
			else{
				//if this is the last path, we then go in opposite direction 
				next_pt_id = (next_pt_id == seg_end) ? seg_st : seg_end;
				break;
			}

			//path_plotter(inward, 0, 0, 0, 1, 0, 0);

			//int lol = 0;


		}
		else{
			//here we are still far away from the re-reouting point
			//so we add next_pt and increment the next_pt_id

			inward.push_back(next_pt);
			points_take_from_that_path++;
			if (next_pt_id == spiral[current_path_id].size() - 1){

				if (current_path_id == spiral.size() - 1){
					break;
				}

				current_path_id++;
				ClipperLib::Path next_path = spiral[current_path_id];
				next_path.swap(current_path);
				next_pt_id = 0;
				//points_take_from_that_path = 0;
			}
			else{
				next_pt_id = next_pt_id + 1;
			}
		}
	}

	VIZ::path_plotter(inward, 0, 0, 1, 1, 0, 0);

	//follow the same logic as for the outward link, but revserse the order 
	//i.e., going in the opposite direction 
	if (reroute_pts.size() > 0){
		contour_end_pt = inward.at(reroute_pts.back());
		reroute_pts.pop_back();
	}
	else{
		//if there is no enough rerouting points i.e., the inward is just
		//one path, then the outward is one path also. the end point is then the p_out
		contour_end_pt = p_out;
	}
	points_take_from_that_path = 0;
	while (true){
		next_pt = current_path[next_pt_id];
		points_take_from_that_path++;
		//check the distance to the before(contour_end_pt)
		ClipperLib::cInt d = dist(ClipperLib::cInt(contour_end_pt.X),
			ClipperLib::cInt(contour_end_pt.Y),
			ClipperLib::cInt(0),
			ClipperLib::cInt(next_pt.X),
			ClipperLib::cInt(next_pt.Y),
			ClipperLib::cInt(0));
		if (d < delta && points_take_from_that_path >current_path.size() / 2){
			outward.push_back(next_pt);

			//decrement the path 
			if (current_path_id > 1){
				//we don't wanna reach path with id=0 because we know
				//that this one is occupied by the inward link
				current_path_id--;
				ClipperLib::Path next_path = spiral[current_path_id];
				next_path.swap(current_path);
				points_take_from_that_path = 0;
			}
			else{
				//add the exit point and exit 
				outward.push_back(contour_end_pt);
				break;
			}


			//reroute by projecttion on current_path 
			int seg_st, seg_end;
			ClipperLib::IntPoint next_proj;
			if (!point_contour_projection(next_pt, current_path, seg_st, seg_end, next_proj, true)){
				PRINT_ERROR("get_fermat():: Can not find the right direction(1)  ");
			}
			else{
				outward.push_back(next_proj);
			}

			int dd_end = dist(next_proj.X, next_proj.Y, ClipperLib::cInt(0),
				current_path[seg_end].X, current_path[seg_end].Y, ClipperLib::cInt(0));
			int dd_st = dist(next_proj.X, next_proj.Y, ClipperLib::cInt(0),
				current_path[seg_st].X, current_path[seg_st].Y, ClipperLib::cInt(0));

			if (dd_end < INT_TOL){
				next_pt_id = seg_end;
			}
			else if (dd_st < INT_TOL){
				next_pt_id = seg_st;
			}
			else{

				ClipperLib::cInt sgn = (next_pt.X - contour_end_pt.X) * (next_proj.Y - contour_end_pt.Y) -
					(next_pt.Y - contour_end_pt.Y) * (next_proj.X - contour_end_pt.X);

				ClipperLib::cInt sgn_st = (next_pt.X - current_path[seg_st].X) * (next_proj.Y - current_path[seg_st].Y) -
					(next_pt.Y - current_path[seg_st].Y) * (next_proj.X - current_path[seg_st].X);

				ClipperLib::cInt sgn_end = (next_pt.X - current_path[seg_end].X) * (next_proj.Y - current_path[seg_end].Y) -
					(next_pt.Y - current_path[seg_end].Y) * (next_proj.X - current_path[seg_end].X);
				//update next_pt_id with the one on the same side as contour_end_pt 
				if (sgn*sgn_st > 0 && sgn*sgn_end < 0){
					next_pt_id = seg_st;
				}
				else if (sgn*sgn_end > 0 && sgn*sgn_st < 0){
					next_pt_id = seg_end;
				}
				else{
					PRINT_ERROR("get_fermat():: Can not find the right direction(2)");
					//TODO: branch out i.e., walk to the right and left and do the same check again 
					//this check fails because the segment could be too small so branching 
					//out solves this
				}

			}


			//update the contour_end_pt
			if (reroute_pts.size() == 0){
				contour_end_pt = p_out;
			}
			else{
				contour_end_pt = inward.at(reroute_pts.back());
				reroute_pts.pop_back();
			}

		}
		else{
			//we are still far from the re-routing point, so we add the next_pt
			//and decrement the point id (and the path if necessary)
			outward.push_back(next_pt);
			points_take_from_that_path++;
			if (next_pt_id == 0){
				if (current_path_id == 0){ break; }
				current_path_id--;
				ClipperLib::Path next_path = spiral[current_path_id];
				next_path.swap(current_path);
				next_pt_id = current_path.size() - 1;
				//points_take_from_that_path = 0;
			}
			else{
				next_pt_id--;
			}
		}

		if (false){
			VIZ::path_plotter(outward, 1, 0, 0, 1, 0, 0);
		}

	}

	VIZ::path_plotter(outward, 1, 0, 0, 1, 0, 0);

	//dot_plotter(p_in.X, p_in.Y, 0, 0, 1, 1);
	//dot_plotter(p_out.X, p_out.Y, 1, 0, 0, 1);
	ClipperLib::Paths fermat;
	fermat.reserve(2);
	fermat.push_back(inward);
	fermat.push_back(outward);

	return fermat;
}

std::vector<ClipperLib::Paths> Spirals::from_contour_to_fermat(
	const std::vector<ClipperLib::Paths>&p, const int delta){

	//take a set of parallel iso-countour and convert it to fermat spiral 
	//only take spirallable regions i.e., p is devided into spirallabel regions 

	//1) create spirals 
	//ClipperLib::Paths spiral = get_spiral(p, delta);

	std::vector<ClipperLib::Paths> spiral = get_spiral(p, delta);

	//visualize the spirals
	//for (int i = 0; i < spiral.size(); i++){
	//	viz_concat_path(spiral[i], false, false, false, true);
	//}

	//2) get fermat spirals out of the created spirals 
	std::vector<ClipperLib::Paths> fermat(spiral.size());
	for (int i = 0; i < spiral.size(); i++){
		fermat[i] = get_fermat(spiral[i], delta);
	}

	//visualize the fermat spirals 
	for (int i = 0; i < fermat.size(); i++){
		double r0, g0, b0, r1, g1, b1;
		if (i == 0){
			r0 = 1.0;
			g0 = 0.0;
			b0 = 1.0;
			r1 = 0.33;
			g1 = 0.67;
			b1 = 1.0;
		}
		else if (i == 1){
			r0 = 0.5;
			g0 = 0.5;
			b0 = 0.0;
			r1 = 0.0;
			g1 = 1.0;
			b1 = 0.0;
		}
		else if (i == 2){
			r0 = 1.0;
			g0 = 0.17;
			b0 = 0.0;
			r1 = 1.0;
			g1 = 0.67;
			b1 = 0.0;
		}
		else if (i == 3){
			r0 = 0.25;
			g0 = 0.25;
			b0 = 0.25;
			r1 = 0.67;
			g1 = 0.67;
			b1 = 0.67;
		}
		else{
			r0 = float(rand()) / float(RAND_MAX);
			g0 = float(rand()) / float(RAND_MAX);
			b0 = float(rand()) / float(RAND_MAX);
			r1 = float(rand()) / float(RAND_MAX);
			g1 = float(rand()) / float(RAND_MAX);
			b1 = float(rand()) / float(RAND_MAX);


		}
		VIZ::path_plotter(fermat[i].at(0), r0, g0, b0, 1, 0, 0, 0);
		VIZ::path_plotter(fermat[i].at(1), r1, g1, b1, 1, 0, 0, 0);

		VIZ::dot_plotter(fermat[i].at(0).front().X, fermat[i].at(0).front().Y, r0, g0, b0, 1);
		VIZ::dot_plotter(fermat[i].at(1).back().X, fermat[i].at(1).back().Y, r1, g1, b1, 1);

	}

	return fermat;


}