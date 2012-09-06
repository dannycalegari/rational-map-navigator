/* magnify.cc magnifying glass */

void rational_map::magnify(){		// magnifying glass
	point p;
	bool mouse_clicked;
	
	p.x=320;
	p.y=320;
	G.magnify_location=p;
	mouse_clicked=false;
	draw_state();
	XFlush(display);
	while(mouse_clicked==false){
		XNextEvent(display, &report);
		switch(report.type) {
			case ButtonPress:
				mouse_clicked=true;
				break;
			case 6:	// MouseMotion?
				p=mouse_location();
				G.magnify_location=p;
				draw_state();
				XFlush(display);
			default:
				break;
		};

	};
};
