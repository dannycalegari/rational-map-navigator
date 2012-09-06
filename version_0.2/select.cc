/* select.cc */

void rational_map::select_and_adjust(){		// select and adjust location of Z/P/V
	point p;
	int mouse_clicks;
	cpx u;
	cpx adjust;
	cvec V;
	
	G.select_type='Z';		// initial values
	G.select_index=-1;
	G.select_location=0.0;
	mouse_clicks=0;
	
	draw_state();
	XFlush(display);
	
	while(G.select_type!='N'){
		XNextEvent(display, &report);
		switch(report.type) {
			case KeyPress:
				if(XLookupKeysym(&report.xkey, 0) == XK_Escape){
					G.select_type='N';
				};
				if(XLookupKeysym(&report.xkey, 0) == XK_s){
					if(G.select_type=='Z'){
						G.select_type='P';
					} else if(G.select_type=='P'){
						G.select_type='V';
					} else if(G.select_type=='V'){
						G.select_type='Z';
					};
				};
				adjust=0.0;
				if(XLookupKeysym(&report.xkey, 0) == XK_Up){
					adjust=0.1*I;
				};
				if(XLookupKeysym(&report.xkey, 0) == XK_Down){
					adjust=-0.1*I;
				};
				if(XLookupKeysym(&report.xkey, 0) == XK_Right){
					adjust=0.1;
				};
				if(XLookupKeysym(&report.xkey, 0) == XK_Left){
					adjust=-0.1;
				};
				if(adjust!=0.0){
					if(G.select_type=='Z'){
						ZERO[G.select_index]=ZERO[G.select_index]+adjust;
						G.select_location=ZERO[G.select_index];
						compute_P_and_Q();
						compute_critical_points(0.000000000000000001);
						compute_critical_values();
					} else if(G.select_type=='P'){
						POLE[G.select_index]=POLE[G.select_index]+adjust;
						G.select_location=POLE[G.select_index];
						compute_P_and_Q();
						compute_critical_points(0.000000000000000001);
						compute_critical_values();
					} else if(G.select_type=='V'){
						V=VALS;
						V[G.select_index]=V[G.select_index]+adjust;
						G.select_location=V[G.select_index];
						flow_VALS_to(V,0.00000001);
					};				
				};
				draw_state();
				XFlush(display);
				break;
			case ButtonPress:
				if(G.select_index==-1){
					p=mouse_location();
					u=point_to_cpx(p);
					switch(G.select_type){
						case 'Z':
							G.select_index=closest_match(u, ZERO);
							if(G.select_index==0){	// can't adjust ZERO[0] which is fixed at 0.
								G.select_index=-1;
							} else {
								G.select_location=ZERO[G.select_index];
							};
							break;
						case 'P':
							G.select_index=closest_match(u, POLE);
							G.select_location=POLE[G.select_index];
							break;
						case 'V':
							p.x=p.x-620;
							u=point_to_cpx(p);
							G.select_index=closest_match(u, VALS);
							G.select_location=VALS[G.select_index];
							break;
						default:
							break;
					};
				} else {
					G.select_index=-1;
				};
				draw_state();
				XFlush(display);
				break;
			default:
				break;
		};
	};
};