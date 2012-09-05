/* points.cc definitions and basic functions of points	*/

struct point{	//
	int x,y;
};

cpx stereo_point(cpx z){
    cpx w;
    double d;
    d=abs(z);
    w.real()=z.real()*3.0/(1.0+d);
    w.imag()=z.imag()*3.0/(1.0+d);
    return(w);
};

cpx inverse_stereo(cpx z){
    cpx w;
    double d;
    d=abs(z);
    w.real()=z.real()/(3.0-d);
    w.imag()=z.imag()/(3.0-d);
    return(w);
};

point cpx_to_point(cpx z){
	cpx w;
	w=stereo_point(z);
    point p;
    p.x = (int) 320.0+100.0*w.real();
    p.y = (int) 320.0-100.0*w.imag();
    return(p);
};

cpx point_to_cpx(point p){
    cpx z;
    z.real()=((double) p.x-320)/100.0;
    z.imag()=((double) 320-p.y)/100.0;
    return(inverse_stereo(z));
};