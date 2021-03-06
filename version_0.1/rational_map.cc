/* rational_map.cc */


class rational_map{
	public:	
	
		// DATA
		vector<complex<double> > Zeros;
		vector<complex<double> > Poles;
		complex<double> M;					// multiplier
		polynomial P,Q;					// R(z) = P(z)/Q(z)
		vector<complex<double> > C;		// critical points
		vector<complex<double> > V;		// critical values
		char ZP;		// is 'Z' for zero or 'P' for pole
		int ZP_index;	// index of selected Zero or Pole
		int V_index;	// index of selected critical Value
		char VF;		// is 'D' for derivative, 'N' for nonlinearity, 'S' for Schwarzian, and 'X' for none
		bool integral_curves;	//
		vector<vector<complex<double> > > PERTURB;		// perturbation matrix; dV/d(Z,P,m)
		vector<complex<double> > STEER;				// steer vector; which way should we move V?
		vector<complex<double> > ADJUST;				// adjust vector; PERTURB(ADJUST) = STEER
		vector<complex<double> > TARGET;			// target vector; where should we move V to?
		vector<transposition> MONODROMY;			// vector of transpositions giving monodromy
		
		// FUNCTIONS
		void initialize(int);				// initialize to "default" values; degree is specified
		void read_map(ifstream &);			// initialize by reading values from a file
		
		complex<double> EVAL(complex<double> );		// evaluate function on complex number
		complex<double> PREIMAGE(complex<double>, complex<double> );	// find inverse near prescribed seed value
		
		void compute_zeros_and_poles();		// determine Zeros and Poles from P and Q
		void compute_coefficients();		// compute coefficients of P and Q
		void compute_C_and_V();				// compute critical points/values
		void adjust_C_and_V();				// adjust values when coefficients move slightly
		rational_map R();					// the function itself!
		rational_map D();					// derivative
		rational_map N();					// nonlinearity
		rational_map Sch();					// Schwarzian
		
		void select_ZP(complex<double>);	// select closest zero/pole to given point
		void adjust_ZP(complex<double>);	// move selected zero/pole by amount
		void select_V(complex<double>);		// select closest critical value to given point
		int select_C(complex<double>);		// index of closest critical point to given point
		void draw_PZCV();					// graphical output routine
		void output_data();					// output data to cout
		
		void compute_monodromy();			// compute monodromy around critical values along ``standard contours''

		void initialize_perturbation_matrix();
		void initialize_target();		// 
		void compute_secant();			// original version of compute_Jacobian 
		void compute_Jacobian();		// computes matrix PERTURB:=dV/d(m,Z,P) 
		void compute_adjust_vector();	// solves PERTURB(ADJUST)=STEER for ADJUST

		void steer_to_target();			// adjust V in a straight line in the direction of TARGET
		void set_target_to_roots_of_unity();
		void set_target_radius(int, double);
		void set_target_argument(int, double);
		void braid_critical_values(int, int, bool);	// braid v_i past [. . . v_{i+/-j}] over (if bool is true) or under (if bool is false) and switch labels
		
		void Mobius();					// adjust R,P by a Mobius transformation
};

void rational_map::initialize(int d){
	int i;
	vector<complex<double> > roots;
	complex<double> w,eta;
	
	w.real()=0;
	w.imag()=TWOPI/(2.0*(double) d);
	eta=exp(w);	// 2dth root of unity
	roots.resize(0);	// initializing roots
	for(i=0;i<d;i++){
		roots.push_back(eta^(2*i));
	};
	roots[0]=5.0;
	Zeros=roots;
	roots.resize(0);
	for(i=0;i<d;i++){
		roots.push_back(eta^(2*i+1));
	};
	roots[0]=-5.0;
	Poles=roots;
	M=1.0;

	compute_coefficients();
	compute_C_and_V();
	initialize_perturbation_matrix();
	initialize_target();
	ZP='Z';			// select zero
	ZP_index=0;		// select zero 0
	V_index=0;		// select critical value 0
	
	VF='X';	// don't draw vector field
	integral_curves=false;	// don't draw integral curves
};

void rational_map::read_map(ifstream &input_file){
	int d,i;
	complex<double> w;
	double r,s;
	input_file >> d;	// read degree
	Zeros.resize(0);
	Poles.resize(0);
	C.resize(0);
	V.resize(0);
	for(i=0;i<d;i++){
		input_file >> r;
		input_file >> s;
		w.real()=r;
		w.imag()=s;
		Zeros.push_back(w);
	};
	for(i=0;i<d;i++){
		input_file >> r;
		input_file >> s;
		w.real()=r;
		w.imag()=s;
		Poles.push_back(w);
	};
	input_file >> r;
	input_file >> s;
	w.real()=r;
	w.imag()=s;
	M=w;
	for(i=0;i<2*d-2;i++){
		input_file >> r;
		input_file >> s;
		w.real()=r;
		w.imag()=s;
		C.push_back(w);
	};
	for(i=0;i<2*d-2;i++){
		input_file >> r;
		input_file >> s;
		w.real()=r;
		w.imag()=s;
		V.push_back(w);
	};
	compute_coefficients();
	adjust_C_and_V();
	initialize_perturbation_matrix();
	initialize_target();
	ZP='Z';			// select zero
	ZP_index=0;		// select zero 0
	V_index=0;		// select critical value 0
	
	VF='X';	// don't draw vector field
	integral_curves=false;	// don't draw integral curves
};

complex<double> rational_map::EVAL(complex<double> z){	// evaluate z
	return(P(z)/Q(z));
};

complex<double> rational_map::PREIMAGE(complex<double> z, complex<double> seed){	// find inverse near prescribed seed value
	complex<double> w,u,v;
	v=seed;
	w=EVAL(v);
	int i;
	i=0;
	while(abs(z-w)>0.00000000001){	// WARNING: accuracy hardcoded! should be able to specify this
		u=D().EVAL(v);		
		v=v-((w-z)/u);	// Newton's method
		w=EVAL(v);
		i++;
		if(i>100){		// if it hasn't found a root quickly, pick new initial value
			i=0;
			real(v)=(double) (100.0*rand() / RAND_MAX)-50.0;
			imag(v)=(double) (100.0*rand() / RAND_MAX)-50.0;
			cout << "root-finding trouble. skipping to new sheet. \n";
		};
	};
	return(v);
};

void rational_map::compute_zeros_and_poles(){
	P.compute_roots();	// P.r is vector of zeros
	Q.compute_roots();	// Q.r is vector of poles
};

void rational_map::compute_coefficients(){
	P=make_polynomial(Zeros,M);
	Q=make_polynomial(Poles,1.0);
};

void rational_map::compute_C_and_V(){	// compute for the first time
	polynomial W;
	int i;
	vector<complex<double> > zv;
	W=Wronskian(P,Q);
	zv=J_T_find_roots(W);
	W.r=zv;
//	W=Wronskian(P,Q);
//	W.compute_roots();	// these are the critical points of P/Q
	C.resize(0);
	V.resize(0);
	for(i=0;i<(int) W.r.size();i++){
		C.push_back(W.r[i]);
		V.push_back(EVAL(W.r[i]));
	};
};

void rational_map::adjust_C_and_V(){	// adjust values, tracking critical values
	polynomial W;
	int i,j;
	complex<double> z;
	vector<complex<double> > zv;
	W=Wronskian(P,Q);
	zv=J_T_find_roots(W);
	W.r=zv;
//	W.compute_roots();	// these are the critical points of P/Q
//	W.compute_roots_with_seed(C);
	vector<complex<double> > L;
	L.resize(0);
	for(i=0;i<(int) W.r.size();i++){
		L.push_back(EVAL(W.r[i]));
	};
	for(i=0;i<(int) C.size();i++){
		j=closest_entry(L,V[i]);	// closest_entry is defined in linear.cc
		C[i]=W.r[j];
		V[i]=EVAL(C[i]);
		W.r.erase(W.r.begin()+j);	// remove this from list of roots
		L.erase(L.begin()+j);

		/* WARNING: finding closest root numerically can (and does) easily lead to
		errors when critical points are too close. It is worth checking the monodromy
		periodically to make sure critical points have not jumped. */
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

rational_map rational_map::R(){		// the function itself
	rational_map R;
	R.P = P;
	R.Q = Q;
	return(R);
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
	for(i=0;i<(int) Zeros.size();i++){
		if(norm(z-Zeros[i])<dist){
			ZP_index=i;
			dist=norm(z-Zeros[i]);
		};
	};
	for(i=0;i<(int) Poles.size();i++){
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
//	compute_C_and_V();
	adjust_C_and_V();
};

void rational_map::select_V(complex<double> z){	// finds critical value closest to z, and sets the value of V_index
	complex<double> w;
	double dist;
	int i;
	V_index=0;
	dist=norm(z-V[0]);
	for(i=0;i<(int) V.size();i++){
		if(norm(z-V[i])<dist){
			V_index=i;
			dist=norm(z-V[i]);
		};
	};
};

int rational_map::select_C(complex<double> z){	// finds critical value closest to z, and sets the value of V_index
	complex<double> w;
	double dist;
	int i;
	int C_index=0;
	dist=norm(z-C[0]);
	for(i=0;i<(int) C.size();i++){
		if(norm(z-C[i])<dist){
			C_index=i;
			dist=norm(z-C[i]);
		};
	};
	return(C_index);
};

void rational_map::output_data(){					// output data to cout
	int i;
	cout << "rational map has the form m.(z-z_0)(z-z_1)...(z-z_{d-1})/(z-p_0)...(z-p_{d-1})\n";
	for(i=0;i<(int) Zeros.size();i++){
		cout << "zero z_" << i << " is " << Zeros[i].real() << " + " << Zeros[i].imag() << " i\n";
	};
	for(i=0;i<(int) Poles.size();i++){
		cout << "pole p_" << i << " is " << Poles[i].real() << " + " << Poles[i].imag() << " i\n";
	};
	cout << "multiplier m is " << M.real() << " + " << M.imag() << "\n";
	for(i=0;i<(int) C.size();i++){
		cout << "critical point " << i << " is " << C[i].real() << " + " << C[i].imag() << " i\n";
	};
	for(i=0;i<(int) V.size();i++){
		cout << "critical value " << i << " is " << V[i].real() << " + " << V[i].imag() << " i\n";
	};
	for(i=0;i<(int) MONODROMY.size();i++){
		cout << "monodromy around critical value " << i << " is " << MONODROMY[i].i << "<->" << MONODROMY[i].j << "\n";
	};
};

void rational_map::compute_monodromy(){			// compute monodromy around critical values along ``standard contours''
	/* Function computes monodromy by steering in preimage space from opposite sides of a
	critical value to zero, computing preimage using Newton's method, seeded to try to be as
	continuous as possible. More robust would be to start at a (preimage of) 0 and steer along
	a ray towards a critical value, and see if we hit a critical point. Note that for generic
	coefficients, all critical points are *simple*, so the monodromy acts as an elementary
	permutation; thus for each critical value, we just need to work out which pair of zeros
	steer into the corresponding critical point.	*/

	int i,j,k;
	complex<double> v,w,d;
	complex<double> I (0.0,1.0);
	point p;
	transposition tau;
		
	MONODROMY.resize(0);
	for(i=0;i<(int) C.size();i++){
		tau.i=-1;						// dummy values
		tau.j=-1;
		MONODROMY.push_back(tau);
	};

	cout << "computing monodromy.\n";
	/* Should check to make sure the argument of R(infty) is not of the form n2\pi/d. 
	If it is, should steer along path slightly kinked away from zero. */
	
	for(i=0;i<(int) V.size();i++){	// ith critical value
		for(j=0;j<(int) Zeros.size();j++){
			w=Zeros[j];
			v=V[i]*(1.0-0.1*I);
			while(abs(EVAL(w)-v)>0.0000000001){	// want to steer w towards v
				d=(D()).EVAL(w);	// derivative at w
				w=w+(0.01*(v-EVAL(w))/d);
				p=complex_to_point(stereo_point(w));
				draw_point(p.x,p.y,0x100410);
				p=complex_to_point(stereo_point(EVAL(w)));
				p.x=p.x+640;
				draw_point(p.x,p.y,0x100410);
			};
			v=V[i]*(1.1);
			while(abs(EVAL(w)-v)>0.0000000001){	// want to steer w towards v
				d=(D()).EVAL(w);	// derivative at w
				w=w+(0.01*(v-EVAL(w))/d);
				p=complex_to_point(stereo_point(w));
				draw_point(p.x,p.y,0x100410);
				p=complex_to_point(stereo_point(EVAL(w)));
				p.x=p.x+640;
				draw_point(p.x,p.y,0x100410);
			};			
			v=V[i]*(1.0+0.1*I);
			while(abs(EVAL(w)-v)>0.0000000001){	// want to steer w towards v
				d=(D()).EVAL(w);	// derivative at w
				w=w+(0.01*(v-EVAL(w))/d);
				p=complex_to_point(stereo_point(w));
				draw_point(p.x,p.y,0x100410);
				p=complex_to_point(stereo_point(EVAL(w)));
				p.x=p.x+640;
				draw_point(p.x,p.y,0x100410);
			};
			v=0.0;
			while(abs(EVAL(w)-v)>0.0000000001){	// want to steer w towards v
				d=(D()).EVAL(w);	// derivative at w
				w=w+(0.01*(v-EVAL(w))/d);
				p=complex_to_point(stereo_point(w));
				draw_point(p.x,p.y,0x100410);
				p=complex_to_point(stereo_point(EVAL(w)));
				p.x=p.x+640;
				draw_point(p.x,p.y,0x100410);
			};

			select_ZP(w); // Zeros[k] is closest zero to w.
			k=ZP_index;
			if(k!=j && ZP=='Z'){
				MONODROMY[i].i=j;
				MONODROMY[i].j=k;
				j=(int) Zeros.size();
			};
		};
	};

	for(i=0;i<(int) MONODROMY.size();i++){
		enforce_order(MONODROMY[i]);
	};
	write_transposition_sequence(MONODROMY);
};

void rational_map::initialize_perturbation_matrix(){
	int i,j;
	vector<complex<double> > COL;
	for(i=0;i<(int) V.size();i++){
		COL.push_back(0.0);
	};
	for(j=0;j<2*(int) Zeros.size()+1;j++){
		PERTURB.push_back(COL);
	};
};

void rational_map::initialize_target(){
	int i;
	TARGET.resize(0);
	for(i=0;i<(int) V.size();i++){
		TARGET.push_back(V[i]);
	};
};


void rational_map::compute_adjust_vector(){ // solves PERTURB(ADJUST)=STEER for ADJUST
	ADJUST=invert_matrix(PERTURB, STEER);
};	// how should we perturb Z, P, m to make V move in the direction STEER?


void rational_map::set_target_to_roots_of_unity(){	// useful to put V in a standard location
	complex<double> w, eta;
	int i;
	
	TARGET.resize(0);
	w.real()=0;
	w.imag()=TWOPI/(double) V.size();
	eta=exp(w);	// 2d-2th root of unity
	for(i=0;i<(int) V.size();i++){
		TARGET.push_back(eta^i);
	};
};

void rational_map::set_target_radius(int i, double t){	// adjusts TARGET[i] to t*TARGET[i]/abs(TARGET[i])
	TARGET[i]=t*TARGET[i]/abs(TARGET[i]);
};

void rational_map::set_target_argument(int i, double t){ // adjusts TARGET[i] to e^{it}abs(TARGET[i])
	complex<double> I;
	I.real()=0.0;
	I.imag()=1.0;
	TARGET[i]=abs(TARGET[i])*exp(t*I);
};

void rational_map::compute_Jacobian(){
	/* For each zero z_i we form a new 2-variable function R(z,l)=R(z-l)/(z-z_i) such that when
	l=z_i we get R(z). Basically, we make the constant z_i into a variable l.
	We also form W(z,l) = Wronskian(P(z-l),Q(z-z_i)) which vanishes at l=z_i when W(z) does.
	
	Then R(c+d,z_i+e) = R(c,z_i) + d dR/dz(c,z_i) + e dR/dl(c,z_i) + quadratic terms.
	Also, W(c+d,z_i+e) = W(c,z_i) + d dW/dz(c,z_i) + e dW/dl(c,z_i) + quadratic terms, so
	if we want W(c+d,z_i+e) = 0 + quadratic terms, we must have d = -e (dW/dl)/(dW/dz)[c,z_i].
	Hence the derivative dv/dl(z_i), where v(l)=R(c(l)), and c(l) is the zero of W(*,l) near c, is
	dR/dl - (dW/dl)*(dR/dz)/(dW/dz) evaluated at c,z_i.
	Explicitly, this is PQR'/W'(c-z_i)^2 - R/(c-z_i).
	
	For each pole p_i we form R(z,l) = R(z-p_i)/(z-l) and W(z,l) = Wr(P(z-p_i),Q(z-l)). We get
	a corresponding formula R/(c-p_i) - PQR'/W'(c-p_i)^2.	*/
	
	int i,j;
	complex<double> z,w;
	vector<complex<double> > COL;
	PERTURB.resize(0);
	polynomial W;
	W=Wronskian(P,Q);
	
	COL=V;
	for(i=0;i<(int) V.size();i++){
		COL[i]=V[i]/M;
	};
	PERTURB.push_back(COL);			// derivative of V with respect to M

	for(i=0;i<(int) Zeros.size();i++){	// correct analytic formula for derivative
		for(j=0;j<(int) C.size();j++){
			w=P.EVAL(C[j])*Q.EVAL(C[j])*(D()).EVAL(C[j])/((W.D()).EVAL(C[j])*(C[j]-Zeros[i])*(C[j]-Zeros[i]));
			w=w-EVAL(C[j])/(C[j]-Zeros[i]);
			COL[j]=w;
		};
		PERTURB.push_back(COL);	
	};

	for(i=0;i<(int) Poles.size();i++){	// correct analytic formula for derivative (or do I have a sign error?)
		for(j=0;j<(int) C.size();j++){
			w=P.EVAL(C[j])*Q.EVAL(C[j])*(D()).EVAL(C[j])/((W.D()).EVAL(C[j])*(C[j]-Poles[i])*(C[j]-Poles[i]));
			w=w-EVAL(C[j])/(C[j]-Poles[i]);
			COL[j]=-w;
		};
		PERTURB.push_back(COL);		// derivative of V with respect to Zeros[i]
	};
};

void rational_map::Mobius(){	// adjust Z,P by a Mobius transformation to prevent clustering if possible.
	int i,j,k,l;
	int ii,jj,kk,ll;
	complex<double> w,ww,I,s,z,zz;	
	double t;
	vector<complex<double> > Z_and_P, original_V;
	complex<double> A,B,C,D,a,b,c,d;
	mmatrix MM,N;
	
	original_V=V;
	z=V[0];
	for(i=0;i<(int) V.size();i++){
		if(norm(V[i])<norm(z)){
			z=V[i];
		};
	};

	Z_and_P = Zeros;
	Z_and_P.insert( Z_and_P.end(), Poles.begin(), Poles.end() );	// list of zeros and poles
	
	ww=cross_ratio(Z_and_P[0],Z_and_P[1],Z_and_P[2],Z_and_P[3]);
	ii=0;
	jj=1;
	kk=2;
	ll=3;
	
	for(i=0;i<(int) Z_and_P.size();i++){
		for(j=0;j<(int) Z_and_P.size();j++){
			for(k=0;k<(int) Z_and_P.size();k++){
				for(l=0;l<(int) Z_and_P.size();l++){
					if((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)!=0){	// if all 4 indices are distinct
						w=cross_ratio(Z_and_P[i],Z_and_P[j],Z_and_P[k],Z_and_P[l]);
						if(norm(w)<norm(ww)){
							ii=i;
							jj=j;
							kk=k;
							ll=l;
							ww=w;
						};
					};
				};
			};
		};
	};
	b=Z_and_P[jj];
	c=Z_and_P[kk];
	d=Z_and_P[ll];
	t=sqrt(abs(ww));
	I.real()=0.0;
	I.imag()=1.0;
	s=exp(t*I);

	MM.A=1.0;
	MM.B=-b;
	MM.C=0.0;
	MM.D=1.0;
	b=operate(MM,b);
	c=operate(MM,c);
	d=operate(MM,d);	// M shifts b to 0
	N.A=1.0/c;
	N.B=0.0;
	N.C=0.0;
	N.D=1.0;
	MM=mult(N,MM);
	b=operate(N,b);
	c=operate(N,c);
	d=operate(N,d);
	N.B=0.0;
	N.D=1.0;
	N.C=(s-d)/(d-d*s);
	N.A=N.C+1.0;
	MM=mult(N,MM);
	N.A=1.0;
	N.B=-0.5;
	N.C=0.0;
	N.D=1.0;
	MM=mult(N,MM);
	N.A=5.0;
	N.B=0.0;
	N.C=0.0;
	N.D=1.0;
	MM=mult(N,MM);
	
	// should take log of MM, and adjust by Z -> Z + epsilon*operate(log(M),Z)

	for(i=0;i<(int) Zeros.size();i++){
		Zeros[i]=operate(MM,Zeros[i]);
	};
	for(i=0;i<(int) Poles.size();i++){
		Poles[i]=operate(MM,Poles[i]);
	};
//  how does M adjust?
	compute_coefficients();
	adjust_C_and_V();
	
	zz=V[0];
	for(i=0;i<(int) V.size();i++){
		if(norm(V[i])<norm(zz)){
			zz=V[i];
		};
	};
	M=M*z/zz;
//	V=original_V;
	compute_coefficients();
	V=original_V;
	adjust_C_and_V();


};

void rational_map::compute_secant(){ 	// how perturbing Z, P and m affects V
	/* original version; computed approximation to Jacobian by secants; reasonably accurate,
	but very slow, because computing new values of C and V after adjusting each variable
	depends on finding a new set of roots for the new Wronskian, via Newton's method, which
	doesn't always converge quickly. 	*/

	int i,j;
	complex<double> z,w;
	vector<complex<double> > COL,VV; 	
	
	VV=V;
	COL=VV;
	for(i=0;i<(int) V.size();i++){
		COL[i]=V[i]/M;
	};
	PERTURB[0]=COL;		// derivative of V with respect to M
	
	for(i=0;i<(int) Zeros.size();i++){
		COL=VV;
		Zeros[i]=Zeros[i]+0.000001;
		compute_coefficients();
		adjust_C_and_V();
		for(j=0;j<(int) V.size();j++){
			COL[j]=(V[j]-COL[j])/0.000001;	// approximate derivative
		};
		Zeros[i]=Zeros[i]-0.000001;
//		compute_coefficients();		this is just *stupid*
//		adjust_C_and_V();
		PERTURB[i+1]=COL;	// derivative of V with respect to Z[i]
	};
	for(i=0;i<(int) Poles.size();i++){
		COL=VV;
		Poles[i]=Poles[i]+0.000001;
		compute_coefficients();
		adjust_C_and_V();
		for(j=0;j<(int) V.size();j++){
			COL[j]=(V[j]-COL[j])/0.000001;	// approximate derivative
		};
		Poles[i]=Poles[i]-0.000001;
//		compute_coefficients();
//		adjust_C_and_V();
		PERTURB[i+Zeros.size()+1]=COL;	// derivative of V with respect to P[i]
	};
	compute_coefficients();
	adjust_C_and_V();
};

