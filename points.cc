/* points.cc complex numbers and points */

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

void draw_point(int, int, long);	// function declaration
