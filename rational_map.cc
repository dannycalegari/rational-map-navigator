/* rational_map.cc */

class rational_map{
	public:	
		vector<complex<double> > Zeros;
		vector<complex<double> > Poles;
		complex<double> M;					// multiplier
		polynomial P,Q;					// R(z) = P(z)/Q(z)
		vector<complex<double> > C;		// critical points
		vector<complex<double> > V;		// critical values
		char ZP;		// is 'Z' for zero or 'P' for pole
		int ZP_index;	// index of ZP
		
		complex<double> EVAL(complex<double> );		// evaluate function on complex number
		void compute_zeros_and_poles();		// determine Zeros and Poles from P and Q
		void compute_coefficients();		// compute coefficients of P and Q
		void compute_C_and_V();				// compute critical points/values
		void select_ZP(complex<double>);	// select closest zero/pole to given point
		void adjust_ZP(complex<double>);	// move selected zero/pole by amount
		void draw_PZCV();					// graphical output routine
		void output_data();					// output data to cout
};

complex<double> rational_map::EVAL(complex<double> z){	// evaluate z
	return(P(z)/Q(z));
};

void rational_map::compute_zeros_and_poles(){
	P.compute_roots();	// P.r is vector of zeros
	Q.compute_roots();	// Q.r is vector of poles
};

void rational_map::compute_coefficients(){
	P=make_polynomial(Zeros,M);
	Q=make_polynomial(Poles,1.0);
};

void rational_map::compute_C_and_V(){
	polynomial W;
	int i;
	W=Wronskian(P,Q);
	W.compute_roots();	// these are the critical points of P/Q
	C.clear();
	V.clear();
	for(i=0;i<W.r.size();i++){
		C.push_back(W.r[i]);
		V.push_back(EVAL(W.r[i]));
	};
};

void rational_map::select_ZP(complex<double> z){	// finds zero or pole closest to z, and sets the value of ZP and ZP_index
	complex<double> w;
	double dist;
	int i;
//	cout << "finding closest zero/pole to " << z.real() << " + " << z.imag() << " i\n";
	ZP='Z';
	ZP_index=0;
	dist=norm(z-Zeros[0]);
	for(i=0;i<Zeros.size();i++){
		if(norm(z-Zeros[i])<dist){
			ZP_index=i;
			dist=norm(z-Zeros[i]);
		};
	};
	for(i=0;i<Poles.size();i++){
		if(norm(z-Poles[i])<dist){
			ZP='P';
			ZP_index=i;
			dist=norm(z-Poles[i]);
		};
	};
};

void rational_map::adjust_ZP(complex<double> z){
	complex<double> w;
	if(ZP=='Z'){	// adjust Z[ZP_index]
		w=Zeros[ZP_index];
		w=stereo_point(w);
		w=w+z;
		w=inverse_stereo(w);
		Zeros[ZP_index]=w;
	} else if(ZP=='P'){		// adjust P[ZP_index]
		w=Poles[ZP_index];
		w=stereo_point(w);
		w=w+z;
		w=inverse_stereo(w);
		Poles[ZP_index]=w;
	};
	compute_coefficients();
	compute_C_and_V();
};

void rational_map::output_data(){					// output data to cout
	int i;
	cout << "rational map has the form m.(z-z_0)(z-z_1)...(z-z_{d-1})/(z-p_0)...(z-p_{d-1})\n";
	for(i=0;i<Zeros.size();i++){
		cout << "zero z_" << i << " is " << Zeros[i].real() << " + " << Zeros[i].imag() << " i\n";
	};
	for(i=0;i<Poles.size();i++){
		cout << "pole p_" << i << " is " << Poles[i].real() << " + " << Poles[i].imag() << " i\n";
	};
	cout << "multiplier m is " << M.real() << " + " << M.imag() << "\n";
	for(i=0;i<C.size();i++){
		cout << "critical point " << i << " is " << C[i].real() << " + " << C[i].imag() << " i\n";
	};
	for(i=0;i<V.size();i++){
		cout << "critical value " << i << " is " << V[i].real() << " + " << V[i].imag() << " i\n";
	};
};