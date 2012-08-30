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

complex<double> cross_ratio(complex<double> a, complex<double> b, complex<double> c, complex<double> d){
	return( ((a-b)*(c-d))/((c-b)*(a-d)) );
};

struct mmatrix{	// mobius matrix
	complex<double> A,B,C,D;
};

complex<double> operate(mmatrix M, complex<double> z){
	return((M.A*z + M.B)/(M.C*z + M.D));
};

mmatrix mult(mmatrix M, mmatrix N){
	mmatrix P;
	P.A = (M.A*N.A) + (M.B*N.C);
	P.B = (M.A*N.B) + (M.B*N.D);
	P.C = (M.C*N.A) + (M.D*N.C);
	P.D = (M.C*N.B) + (M.D*N.D);
	return(P);
};