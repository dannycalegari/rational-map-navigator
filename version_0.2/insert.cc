/* insert.cc	routine for inserting zero/pole pair */

void rational_map::insert_zp(){		// insert zero/pole pair
	point p;
	int points_selected;
	
	points_selected=0;
	draw_state();
	XFlush(display);
	while(points_selected<2){
		XNextEvent(display, &report);
		switch(report.type) {
			case ButtonPress:
				p=mouse_location();
				points_selected++;
				if(points_selected==1){
					G.insert_point=point_to_cpx(p);
					G.zero_point=G.insert_point+0.01;
					ZERO.push_back(G.zero_point);
					POLE.push_back(2.0*G.insert_point-G.zero_point);
					CRIT.resize(0);
					VALS.resize(0);
					compute_P_and_Q();		
					compute_critical_points(0.00000001);	
					compute_critical_values();
				};
				draw_state();
				XFlush(display);
				break;
			case 6:	// MouseMotion?
				if(points_selected==1){
					p=mouse_location();
					G.zero_point=point_to_cpx(p);
					ZERO[ZERO.size()-1]=(G.zero_point);
					POLE[POLE.size()-1]=(2.0*G.insert_point-G.zero_point);
					compute_P_and_Q();		
					compute_critical_points(0.00000001);	
					compute_critical_values();
					draw_state();
					XFlush(display);
				};
			default:
				break;
		};
	};
};
