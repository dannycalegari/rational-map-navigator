/*******************************************************
 * graphics.cc - sets up window stuff for main program *
 *                                                     *
 * uses standard Xlib stuff, because I'm too perverse  *
 * to use a GUI toolkit. the main reason is that       *
 * they keep upgrading to newer and newer libraries    *
 * which are not backward compatible, so that my       *
 * programs always start to break. boohoohoo.          *
 *                                                     *
 * Danny Calegari 12/17/2000                           *
 *******************************************************/

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

Display *display;
int screen_num;
unsigned int display_width, display_height;
XEvent report;
GC gc;
Window win;
int border_width = 4;
unsigned int width, height;
XFontStruct * font;

void setup_graphics(void){
	display=XOpenDisplay(NULL);
	display_width = DisplayWidth(display, screen_num);
	display_height = DisplayHeight(display, screen_num);
	screen_num = DefaultScreen(display);  
	width = 1280;
	height = 680;
	win = XCreateSimpleWindow(display, RootWindow(display, screen_num), 0, 0, width, 
		height, border_width, BlackPixel(display, screen_num), WhitePixel(display, screen_num));
	XSelectInput(display, win, ExposureMask | KeyPressMask | ButtonPressMask);
	gc = DefaultGC(display, screen_num);
	XSetForeground(display, gc, BlackPixel(display, screen_num));
	XSetBackground(display, gc, WhitePixel(display, screen_num));
  
	XMapWindow(display, win);
}

void setup_font(void){
    const char * fontname = "-*-georgia-*-r-*-*-14-*-*-*-*-*-*-*";
 //   const char * fontname = "-*-times-*-r-*-*-16-*-*-*-*-*-*-*";

    font = XLoadQueryFont (display, fontname);
    /* If the font could not be loaded, revert to the "fixed" font. */
    if (! font) {
        font = XLoadQueryFont (display, "fixed");
        cout << "couldn't find font!\n";
    }
    XSetFont (display, gc, font->fid);
}


void erase_field(void){
	XClearWindow(display, win);
}



void draw_line(int x1, int y1, int x2, int y2){
        XSetForeground(display, gc, (long) 0);
        XSetLineAttributes(display, gc, 2, LineSolid, 1, 1);
        XDrawLine(display, win, gc, x1, y1, x2, y2);
};

void draw_faint_line(int x1, int y1, int x2, int y2){
        XSetForeground(display, gc, (long) 0xDDDDDD);
        XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
        XDrawLine(display, win, gc, x1, y1, x2, y2);
};

void draw_thin_line(int x1, int y1, int x2, int y2, long col){
        XSetForeground(display, gc, (long) col);
        XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
        XDrawLine(display, win, gc, x1, y1, x2, y2);
};

void draw_circle(int x, int y, int r, long col){
        XSetForeground(display, gc, col);
        XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
        XSetFillStyle(display, gc, FillSolid);
        XDrawArc(display, win, gc, x-r, y-r, 2*r, 2*r, 0, 23040);
};

void draw_thick_circle(int x, int y, int r, long col){
        XSetForeground(display, gc, col);
        XSetLineAttributes(display, gc, 4, LineSolid, 1, 1);
        XDrawArc(display, win, gc, x-r, y-r, 2*r, 2*r, 0, 23040);
};

void draw_faint_circle(int x, int y, int r){
        XSetForeground(display, gc, (long) 0xDDDDDD);
        XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
        XDrawArc(display, win, gc, x-r, y-r, 2*r, 2*r, 0, 23040);
};

void draw_point(int x, int y, long col){
        XSetForeground(display, gc, col);
		XDrawPoint(display, win, gc, x, y);
};

void draw_vector_field(rational_map R, long col){
	int x, y, xx, yy;
	point p;
	complex<double> z,w;
	for(x=20;x<620;x=x+8){
		for(y=20;y<620;y=y+8){
			p.x=x;
			p.y=y;
			z=point_to_complex(p);
			if(abs(z)<3.0){
				z=inverse_stereo(z);
				w=R.EVAL(z);
				xx = (int) (8.0*w.real())/(1.0+abs(w));
				yy = (int) (8.0*w.imag())/(1.0+abs(w));
				draw_thin_line(x,y,x+xx,y+yy, col);
			};
		};
	};
};

void draw_integral_curves(rational_map R, long col){
	int x, y;
	point p,q;
	complex<double> z,w,zz;
	int i;
	for(x=0;x<620;x=x+10){
		for(y=0;y<620;y=y+10){
			p.x=x;
			p.y=y;
			z=point_to_complex(p);
			if(abs(z)<3.0){
				z=inverse_stereo(z);	// OK, got initial, viable point
				for(i=0;i<100;i++){
					w=R.EVAL(z);
					z=z+w/(1.0+abs(w)*100.0);	// step size too small?
		//			cout << z.real() << " + " << z.imag() << "\n";
					zz=stereo_point(z);
					if(abs(zz)<3.0){
						q=complex_to_point(stereo_point(z));
						draw_thin_line(p.x,p.y,q.x,q.y, col);
					};
					p=q;
				};
			};
		};
	};
};

void draw_real_curves(rational_map R, long col){
	int x, y;
	point p;
	complex<double> z,w,zz;
	int i;
	for(x=0;x<620;x=x+1){
		for(y=0;y<620;y=y+1){
			p.x=x;
			p.y=y;
			z=point_to_complex(p);
			if(abs(z)<3.0){
				z=inverse_stereo(z);	// OK, got initial, viable point
				w=R.EVAL(z);
				i=(int) (1000.0+(w.imag()*0.4)) % 2;
				if(i==0){
					draw_point(p.x, p.y, col);
				};
			};
		};
	};
};

void rational_map::draw_PZCV(){	// graphical output routine
	int i;
	point p;
	stringstream T;
	string S;
	complex<double> z;
	draw_line(640,0,640,640);
	draw_line(0,640,1280,640);
	draw_faint_line(310,320,330,320);
	draw_faint_line(320,330,320,310);
	draw_faint_line(950,320,970,320);
	draw_faint_line(960,330,960,310);
	draw_faint_circle(320,320,300);
	draw_faint_circle(960,320,300);
	switch(VF){
		case 'D':
			if(integral_curves){
				draw_integral_curves(D(), 0xBBCCDD); // draw derivative
			} else {
				draw_vector_field(D(), 0x99AABB); // draw derivative
			};
			break;
		case 'N':
			if(integral_curves){
				draw_integral_curves(N(), 0xBBCCDD); // draw nonlinearity
			} else {
				draw_vector_field(N(), 0x99AABB); // draw nonlinearity
			};
			break;
		case 'S':
			if(integral_curves){
				draw_integral_curves(Sch(), 0xBBCCDD); // draw Schwarzian
			} else {
				draw_vector_field(Sch(), 0x99AABB); // draw Schwarzian
			};
		case 'X':
			break;
		default:
			break;
	};

	for(i=0;i<(int) Zeros.size();i++){
		p=complex_to_point(stereo_point(Zeros[i]));
		draw_thick_circle(p.x,p.y,1,(long) 0xFF0000);
	};
	for(i=0;i<(int) Poles.size();i++){
		p=complex_to_point(stereo_point(Poles[i]));
		draw_thick_circle(p.x,p.y,1,(long) 0x0000FF);
	};
	
	if(ZP=='Z'){
		p=complex_to_point(stereo_point(Zeros[ZP_index]));
	} else {
		p=complex_to_point(stereo_point(Poles[ZP_index]));
	};
	draw_circle(p.x,p.y,6,(long) 0x000000);	// big circle around selected Z/P
	
	p=complex_to_point(stereo_point(V[V_index]));
	p.x=p.x+640;
	draw_circle(p.x,p.y,6,(long) 0x000000);	// big circle around selected V
	
	for(i=0;i<(int) C.size();i++){
		p=complex_to_point(stereo_point(C[i]));
		draw_thick_circle(p.x,p.y,1,(long) 0x00FF00);
		T << i;
		S=T.str();
		XDrawString(display,win,gc,p.x+3,p.y-3,S.c_str(),strlen(S.c_str()));
		T.str("");
	};
	for(i=0;i<(int) V.size();i++){
		p=complex_to_point(stereo_point(V[i]));
		p.x=p.x+640;
		draw_thick_circle(p.x,p.y,1,(long) 0x3377AA);
		T << i;
		S=T.str();
		XDrawString(display,win,gc,p.x+3,p.y-3,S.c_str(),strlen(S.c_str()));
		T.str("");
	};
    XSetForeground(display, gc, (long) 0x000000);

	if(ZP=='Z'){
		S="Current selected point is a zero, located at ";
		z=Zeros[ZP_index];
	} else if(ZP=='P') {
		S="Current selected point is a pole, located at ";
		z=Poles[ZP_index];
	};
	T << z.real() << " + " << z.imag() << " i";
	S=S+T.str();
	T.str("");

//	cout << S;
	XDrawString(display,win, gc,120,665,S.c_str(),strlen(S.c_str()));
	switch(VF){
		case 'D':
			S="Drawing derivative.";
			break;
		case 'N':
			S="Drawing nonlinearity.";
			break;
		case 'S':
			S="Drawing Schwarzian.";
			break;
		case 'X':
			S="Not drawing vector field.";
			break;
		default:
			break;	
	};
	if(integral_curves){
		S=S+" Integral curves on.";
	} else {
		S=S+" Integral curves off.";
	};
	XDrawString(display, win, gc, 800,665,S.c_str(),strlen(S.c_str()));
};

point mouse_location(){
//    Bool result;
    Window window_returned;
    int root_x, root_y;
    int win_x, win_y;
    unsigned int mask_return;
    point p;
    
	XQueryPointer(display, win, &window_returned,
                &window_returned, &root_x, &root_y, &win_x, &win_y,
                &mask_return);
    p.x=win_x;
    p.y=win_y;
    return(p);
};

double norm(vector<complex<double> > v){
	int i;
	double t;
	t=0;
	for(i=0;i<(int) v.size();i++){
		t=t+abs(v[i]);
	};
	return(t);
};

void rational_map::steer_to_target(){
	// adjusts zeros/poles to move critical values in a "straight" line to roots of unity
	
	int i,j;
	complex<double> w, eta;
	vector<complex<double> > PROX;
	double SPEED;

	for(i=0;i<(int) V.size();i++){
		PROX.push_back(1.0);
	};
	
	j=0;
	SPEED=0.03;
	STEER=PROX;
	while(norm(PROX)>0.01){
		if(j%4==0){
			erase_field();
			draw_PZCV();
			XFlush(display);
		};
		
		for(i=0;i<(int) V.size();i++){
			PROX[i]=((TARGET[i]-V[i]));
			if(abs(PROX[i])<0.1){	// if close enough,
				STEER[i]=PROX[i]/SPEED;
			} else if(abs(PROX[i])<0.2){
				STEER[i]=PROX[i]*5.0;
			} else	{
				STEER[i]=PROX[i]*(1.0+ 1.0/abs(PROX[i]));
				if(abs(STEER[i])>10.0){	// absolute speed limit
					STEER[i]=STEER[i]*10.0/abs(STEER[i]);
				};
			};
		};		
		compute_perturbation_matrix();
		compute_adjust_vector();
		M=M+SPEED*ADJUST[0];
		for(i=0;i<(int) Zeros.size();i++){
			Zeros[i]=Zeros[i]+SPEED*ADJUST[i+1];
		};
		for(i=0;i<(int) Poles.size();i++){
			Poles[i]=Poles[i]+SPEED*ADJUST[i+Zeros.size()+1];
		};
		compute_coefficients();

		adjust_C_and_V();
		j++;
	};
};

void graphics_routine(rational_map &R, bool &finished){
	// function stereo_point takes complex plane to disk of radius 2, stereographically
	// function inverse_stereo takes disk of radius 2 to complex plane
	complex<double> z;
	point p;
		
	erase_field();
	R.draw_PZCV();

	XFlush(display);
	XNextEvent(display, &report);
	switch (report.type) {
		case ButtonPress:
			p=mouse_location();
			if(p.x<640){
				z=point_to_complex(p);
				R.select_ZP(inverse_stereo(z));
			} else {
				p.x=p.x-640;
				z=point_to_complex(p);
				R.select_V(inverse_stereo(z));
			};
		case KeyPress:
			if(XLookupKeysym(&report.xkey, 0) == XK_Right){		// adjust Z or P
				z.real()=0.04;
				z.imag()=0.00;
				R.adjust_ZP(z);
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_Left){		// adjust Z or P
				z.real()=-0.04;
				z.imag()=0.00;
				R.adjust_ZP(z);
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_Up){		// adjust Z or P
				z.real()=0.00;
				z.imag()=0.04;
				R.adjust_ZP(z);
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_Down){		// adjust Z or P
				z.real()=0.00;
				z.imag()=-0.04;
				R.adjust_ZP(z);				
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_o){		// output data
				R.output_data();
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_v){		// toggle vector field D/N/S
				switch(R.VF){
					case 'X':
						R.VF='D';
						break;
					case 'D':
						R.VF='N';
						break;
					case 'N':
						R.VF='S';
						break;
					case 'S':
						R.VF='X';
						break;
					default:
						R.VF='X';
						break;	
				};
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_i){		// toggle integral curves
				if(R.integral_curves){
					R.integral_curves=false;
				} else {
					R.integral_curves=true;
				};
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_m){		// compute monodromy
				R.compute_monodromy();
				while(1){
					XNextEvent(display, &report);
					if(report.type!=0){	// waiting to interrupt to outer loop
						break;
					};
				};
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_s){		// steer to roots of unity
				R.set_target_to_roots_of_unity();
				R.steer_to_target();
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_a){
				R.set_target_radius(R.V_index,abs(R.V[R.V_index])*1.2);
				R.steer_to_target();
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_z){		
				R.set_target_radius(R.V_index,abs(R.V[R.V_index])/1.2);
				R.steer_to_target();
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_x){
				R.set_target_argument(R.V_index,arg(R.V[R.V_index])+0.1);
				R.steer_to_target();
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_c){
				R.set_target_argument(R.V_index,arg(R.V[R.V_index])-0.1);
				R.steer_to_target();
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_q){		// quit
				finished=true;
				XCloseDisplay(display);
				exit(0);
				break;
			};
			default:
				break;
	};
};

