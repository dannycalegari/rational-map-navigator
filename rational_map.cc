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
		char VF;		// is 'D' for derivative, 'N' for nonlinearity, 'S' for Schwarzian, and 'X' for none
				
		complex<double> EVAL(complex<double> );		// evaluate function on complex number
		void compute_zeros_and_poles();		// determine Zeros and Poles from P and Q
		void compute_coefficients();		// compute coefficients of P and Q
		void compute_C_and_V();				// compute critical points/values
		rational_map D();					// derivative
		rational_map N();					// nonlinearity
		rational_map Sch();					// Schwarzian
		
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

rational_map operator+(rational_map R, rational_map S){
	rational_map T;
	T.P=(R.P*S.Q)+(S.P*R.Q);
	T.Q=R.Q*S.Q;
	return(T);
};

rational_map operator-(rational_map R, rational_map S){
	rational_map T;
	T.P=(R.P*S.Q)-(S.P*R.Q);
	T.Q=R.Q*S.Q;
	return(T);
};

rational_map operator*(rational_map R, rational_map S){
	rational_map T;
	T.P=R.P*S.P;
	T.Q=R.Q*S.Q;
	return(T);
};

rational_map operator/(rational_map R, rational_map S){
	rational_map T;
	T.P=R.P*S.Q;
	T.Q=R.Q*S.P;
	return(T);
};

rational_map operator*(rational_map R, polynomial P){
	rational_map S;
	S.P=R.P*P;
	S.Q=R.Q;
	return(S);
};

rational_map operator*(rational_map R, complex<double> z){
	polynomial P;
	P=monomial(z,0);
	return(R*P);
};

rational_map rational_map::D(){		// derivative
	rational_map R;
	R.P=Wronskian(P,Q);
	R.Q=Q*Q;
	return(R);
};

rational_map D(rational_map R){
	rational_map S;
	S.P=Wronskian(R.P,R.Q);
	S.Q=R.Q*R.Q;
	return(S);
};

rational_map rational_map::N(){		// nonlinearity
	/* Nonlinearity is (log f')' = f''/f'
	 if g=f', nonlinearity is g'/g
	 if g=P/Q, nonlinearity is (P'Q-Q'P)*Q/P*Q^2 = P'/P - Q'/Q = (P'Q-Q'P)/PQ
	 if f=P/Q so g=(P'Q-Q'P)/Q^2 then N(f) = (P''Q-Q''P)/(P'Q-Q'P) - 2Q'/Q
	 = [(P''Q-Q''P)*Q - 2Q'(P'Q-Q'P)] / Q(P'Q-Q'P)
	 = (P''Q^2 - Q''PQ - 2P'Q'Q - 2(Q')^2P) / Q(P'Q-Q'P) */
	rational_map R;
	R.P=Q*(Q*P.D().D() - P*Q.D().D())-Q.D()*Wronskian(P,Q)-Q.D()*Wronskian(P,Q);
	R.P.a.pop_back();	// (last coefficient is zero)
	R.Q=Q*Wronskian(P,Q);
	return(R);
};

rational_map Schwarzian(rational_map R){		// returns the Schwarzian of the rational map
	// Schwarzian is (f''/f')' - (f''/f')^2/2

	return( D( D(D(R))/D(R) ) - ( D(D(R))*D(D(R))*0.5 / (D(R)*D(R)) ) );
};	

rational_map rational_map::Sch(){
	/* S.P = -2P'''PQQ' + 6PQP''Q'' - 6PP''Q'Q' - 3QQP''P'' - 2PQQ'''P' -
			6QP'P'Q'' + 6PP'Q'Q'' + 2P'''QQP' + 6QP'P''Q' - 3PPQ''Q'' + 2PPQ'''Q'
	   S.Q = 2(QP' - PQ')^2
	*/
	rational_map R;
	R.P = ((P.D().D().D())*Q.D()*P*Q)*(-2.0) + ((P.D().D())*(Q.D().D())*P*Q)*6.0 + ((P.D().D())*Q.D()*Q.D()*P)*(-6.0) +
		((P.D().D())*(P.D().D())*Q*Q)*(-3.0) + ((Q.D().D().D())*P.D()*P*Q)*(-2.0) + ((Q.D().D())*P.D()*P.D()*Q)*(-6.0) +
		((Q.D().D())*Q.D()*P.D()*P)*6.0 + ((P.D().D().D())*P.D()*Q*Q)*2.0 + ((P.D().D())*P.D()*Q.D()*Q)*(6.0) +
		((Q.D().D())*(Q.D().D())*P*P)*(-3.0) + ((Q.D().D().D())*Q.D()*P*P)*2.0;
	R.Q = (Wronskian(P,Q)*Wronskian(P,Q))*2.0;
	R.P.a.pop_back();	// (last 4 coefficients are zero, because of cancellation)
	R.P.a.pop_back();	
	R.P.a.pop_back();	
	R.P.a.pop_back();	

//	cout << R.P.a[R.P.a.size()-1].real() << " + " << R.P.a[R.P.a.size()-1].imag() << "\n";
//	cout << R.Q.a[R.Q.a.size()-1].real() << " + " << R.Q.a[R.Q.a.size()-1].imag() << "\n";
	return(R);
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