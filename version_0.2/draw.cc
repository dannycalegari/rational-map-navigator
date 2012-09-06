/* draw.cc draw state */

void rational_map::draw_state(){
	int i;
	cpx v;
	point p,q;
	stringstream T;
	
	erase_field();
	
	// special: help screen
	
	if(G.doing=='H'){
		p.x=30;
		p.y=60;
		
		T << "Welcome to rational map navigator. ";
		T << "This is the help screen.";
		p.y=p.y+30;
		draw_text(p,T,0x000000);
		T.str("");				
		T << "On LHS: zeros are red, poles are blue, critical points are green.";
		p.y=p.y+60;
		draw_text(p,T,0x000000);
		T.str("");	
		T << "On RHS: critical values are purple.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");			
		T << "Type [f] to flow to roots of unity.";
		p.y=p.y+60;
		draw_text(p,T,0x000000);
		T.str("");			
		T << "Type [l] to toggle labels on/off.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");	
		T << "Type [h] to toggle help screen on/off.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");	
		T << "Type [i] to insert a zero/pole pair.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");			
		T << "Type [s] to select/adjust Z/P/V.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");			
		T << "Type [m] for magnifying glass.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");					
		T << "Type [w] to write to file.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");
		T << "Type [q] to quit.";
		p.y=p.y+40;
		draw_text(p,T,0x000000);
		T.str("");	
		XFlush(display);
	} else {
	
	// draw outer circles

		p.x=320;
		p.y=320;
		draw_circle(p,300,0xAAAAAA);
		p.x=p.x+620;
		draw_circle(p,300,0xAAAAAA);
	
	// draw Z/P/C/V
	
		for(i=0;i<(int) ZERO.size();i++){
			p=cpx_to_point(ZERO[i]);
			draw_concentric_circles(p,2,0x550000);
			if(G.labels_on){
				draw_label(p,i,0x550000);
			};
		};
		for(i=0;i<(int) POLE.size();i++){
			p=cpx_to_point(POLE[i]);
			draw_concentric_circles(p,2,0x000055);
			if(G.labels_on){
				draw_label(p,i,0x000055);
			};
		};
		for(i=0;i<(int) CRIT.size();i++){
			p=cpx_to_point(CRIT[i]);
			draw_concentric_circles(p,2,0x005500);
			if(G.labels_on){
				draw_label(p,i,0x005500);
			};
		};
		for(i=0;i<(int) VALS.size();i++){
			p=cpx_to_point(VALS[i]);
			p.x=p.x+620;
			draw_concentric_circles(p,2,0x550055);
			if(G.labels_on){
				draw_label(p,i,0x550055);
			};
		};
		
	// write state
	
		T << "Degree " << deg(P) << ".  ";
		T << "Type [h] for help screen.  ";
		
		if(G.labels_on){
			T << "Labels are on.  ";
		} else {
			T << "Labels are off.  ";
		};
			
		switch(G.doing){
			case 'F':
				T << "Flowing critical values to roots of unity. Distance is " << G.distance << ".  ";
				break;
			case 'U':
				T << "User control.  ";
				break;
			case 'I':
				T << "Insert zero/pole. Use mouse to pick location.  ";
				break;
			case 'M':
				T << "Magnifying glass. Move mouse over location to magnify.  ";
				break;
			case 'S':
				T << "Select Z/P/V. Use mouse to select and move point.  ";
			default:
				break;
		};
		
		p.x=30;
		p.y=660;
		draw_text(p,T,0x000000);
		T.str("");

	// special: magnifying glass mode; erase circle and redraw some points 
	
		if(G.doing=='M'){
			erase_circle(G.magnify_location,196);
			draw_circle(G.magnify_location,196,0x000000);
			for(i=0;i<(int) ZERO.size();i++){
				p=cpx_to_point(ZERO[i]);
				q.x=p.x-G.magnify_location.x;
				q.y=p.y-G.magnify_location.y;
				if(q.x*q.x+q.y*q.y < 784){
					p.x=G.magnify_location.x+q.x*7;
					p.y=G.magnify_location.y+q.y*7;
					draw_concentric_circles(p,2,0x550000);
					if(G.labels_on){
						draw_label(p,i,0x550000);
					};
				};
			};
			
			for(i=0;i<(int) POLE.size();i++){
				p=cpx_to_point(POLE[i]);
				q.x=p.x-G.magnify_location.x;
				q.y=p.y-G.magnify_location.y;
				if(q.x*q.x+q.y*q.y < 784){
					p.x=G.magnify_location.x+q.x*7;
					p.y=G.magnify_location.y+q.y*7;
					draw_concentric_circles(p,2,0x000055);
					if(G.labels_on){
						draw_label(p,i,0x000055);
					};
				};
			};
			for(i=0;i<(int) CRIT.size();i++){
				p=cpx_to_point(CRIT[i]);
				q.x=p.x-G.magnify_location.x;
				q.y=p.y-G.magnify_location.y;
				if(q.x*q.x+q.y*q.y < 784){
					p.x=G.magnify_location.x+q.x*7;
					p.y=G.magnify_location.y+q.y*7;
					draw_concentric_circles(p,2,0x005500);
					if(G.labels_on){
						draw_label(p,i,0x005500);
					};
				};
			};
			for(i=0;i<(int) VALS.size();i++){
				p=cpx_to_point(VALS[i]);
				p.x=p.x+620;
				q.x=p.x-G.magnify_location.x;
				q.y=p.y-G.magnify_location.y;
				if(q.x*q.x+q.y*q.y < 784){
					p.x=G.magnify_location.x+q.x*7;
					p.y=G.magnify_location.y+q.y*7;
					draw_concentric_circles(p,2,0x550055);
					if(G.labels_on){
						draw_label(p,i,0x550055);
					};
				};
			};
		};
	};
	XFlush(display);
};
