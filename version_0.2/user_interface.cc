/* user_interface.cc */

void rational_map::user_interface(){
	point p;
	bool finished;
	cvec V;
	int i;
	cpx u;
	ofstream output_file;
	braid b;
	
	finished=false;
	while(finished==false){
		XNextEvent(display, &report);
		switch(report.type) {
			case ButtonPress:
				p=mouse_location();
				break;
			case KeyPress:
				if(XLookupKeysym(&report.xkey, 0) == XK_f){	// flow to roots of unity
					V.resize(0);
					for(i=0;i<(int) VALS.size();i++){
						u=exp(TWOPI*I*(double) i/(double) VALS.size());
						V.push_back(u);
					};
					G.doing='F';
					flow_VALS_to(V,0.00000001);
					G.doing='U';
					draw_state();
					break;
				};
				if(XLookupKeysym(&report.xkey, 0) == XK_q){ // quit           
                    finished=true;
                    XCloseDisplay(display);
                    exit(0);
                    break;
                };
                if(XLookupKeysym(&report.xkey, 0) == XK_l){	// toggle labels on/off
                	G.labels_on=1-G.labels_on;
                	draw_state();
                	break;
                };
                if(XLookupKeysym(&report.xkey, 0) == XK_h){	// help screen
                	if(G.doing=='U'){
	                	G.doing='H';
	                } else {
	                	G.doing='U';
	                };
                	draw_state();
                	break;
                };
                if(XLookupKeysym(&report.xkey, 0) == XK_i){	// insert zero/pole
                	G.doing='I';
                	insert_zp();
                	G.doing='U';
                	draw_state();
                	break;
                };
                if(XLookupKeysym(&report.xkey, 0) == XK_m){ // magnifying glass
                	G.doing='M';
                	magnify();
                	G.doing='U';
                	draw_state();
                	break;
                };
                if(XLookupKeysym(&report.xkey, 0) == XK_w){ // write to file
                	output_file.open("output_file.txt");
					write_to_file(output_file);
					output_file.close();
					break;
				};
				if(XLookupKeysym(&report.xkey, 0) == XK_b){	// do braid
					b.i=2;
					b.j=-5;
					b.over=false;
					do_braid(b);
					break;
				};
            default:
            	break;
        };
    };
};