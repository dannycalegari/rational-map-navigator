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
    const char * fontname = "-*-times-*-r-*-*-14-*-*-*-*-*-*-*";
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

struct point{
	int x,y;
};

point complex_to_point(complex<double> z){
	point p;
	p.x = (int) 320.0+100.0*z.real();
	p.y = (int) 320.0-100.0*z.imag();
	return(p);
};

complex<double> point_to_complex(point p){
	complex<double> z;
	z.real()=((double) p.x-320)/100.0;
	z.imag()=((double) 320-p.y)/100.0;
	return(z);
};

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

void draw_circle(int x, int y, int r, long col){
        XSetForeground(display, gc, col);
        XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
        XDrawArc(display, win, gc, x-r, y-r, 2*r, 2*r, 0, 23040);
};

void draw_faint_circle(int x, int y, int r){
        XSetForeground(display, gc, (long) 0xDDDDDD);
        XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
        XDrawArc(display, win, gc, x-r, y-r, 2*r, 2*r, 0, 23040);
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
	draw_faint_line(970,320,990,320);
	draw_faint_line(980,330,980,310);
	draw_faint_circle(320,320,300);
	draw_faint_circle(980,320,300);
	for(i=0;i<Zeros.size();i++){
		p=complex_to_point(stereo_point(Zeros[i]));
		draw_circle(p.x,p.y,1,(long) 0xFF0000);
	};
	for(i=0;i<Poles.size();i++){
		p=complex_to_point(stereo_point(Poles[i]));
		draw_circle(p.x,p.y,1,(long) 0x0000FF);
	};
	if(ZP=='Z'){
		p=complex_to_point(stereo_point(Zeros[ZP_index]));
	} else {
		p=complex_to_point(stereo_point(Poles[ZP_index]));
	};
	draw_circle(p.x,p.y,3,(long) 0x000000);	// big circle around selected Z/P
	
	for(i=0;i<C.size();i++){
		p=complex_to_point(stereo_point(C[i]));
		draw_circle(p.x,p.y,1,(long) 0x00FF00);
	};
	for(i=0;i<V.size();i++){
		p=complex_to_point(stereo_point(V[i]));
		p.x=p.x+640;
		draw_circle(p.x,p.y,1,(long) 0x3377AA);
	};
	if(ZP=='Z'){
		S="Current selected point is a zero, located at ";
		z=Zeros[ZP_index];
	} else if(ZP=='P') {
		S="Current selected point is a pole, located at ";
		z=Poles[ZP_index];
	};
	T << z.real() << " + " << z.imag() << " i";
	S=S+T.str();
//	cout << S;
    XSetForeground(display, gc, (long) 0x000000);
	XDrawString(display,win, gc,20,665,S.c_str(),strlen(S.c_str()));
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


void graphics_routine(rational_map &R, bool &finished){
	// function stereo_point takes complex plane to disk of radius 2, stereographically
	// function inverse_stereo takes disk of radius 2 to complex plane
	complex<double> z;
	point p;
	bool result;
		
	erase_field();
	R.draw_PZCV();

	XFlush(display);
	XNextEvent(display, &report);
	switch (report.type) {
		case ButtonPress:
			p=mouse_location();
			z=point_to_complex(p);
			R.select_ZP(inverse_stereo(z));
		case KeyPress:
			if(XLookupKeysym(&report.xkey, 0) == XK_Right){		// adjust Z or P
				z.real()=0.02;
				z.imag()=0.00;
				R.adjust_ZP(z);
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_Left){		// adjust Z or P
				z.real()=-0.02;
				z.imag()=0.00;
				R.adjust_ZP(z);
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_Up){		// adjust Z or P
				z.real()=0.00;
				z.imag()=0.02;
				R.adjust_ZP(z);
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_Down){		// adjust Z or P
				z.real()=0.00;
				z.imag()=-0.02;
				R.adjust_ZP(z);				
				break;
			};
			if(XLookupKeysym(&report.xkey, 0) == XK_o){		// output data
				R.output_data();
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