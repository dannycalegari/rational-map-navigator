/* roots.cc function to compute roots of complex polynomial */

struct real_polynomial{
	vector<double> a;
};

real_polynomial D(real_polynomial P){
	real_polynomial Q;
	Q.a.resize(0);
	int i;
	for(i=1;i<(int) P.a.size();i++){
	 	Q.a.push_back(P.a[i]*(double) i);
	};
	return(Q);
};

double EVAL(real_polynomial P, double t){
	int i;	
	double s;
	s=0.0;
	for(i=(int) P.a.size()-1;i>=0;i--){
		s=(s*t) + P.a[i];
	};
	return(s);	
};

double find_real_positive_simple_root(real_polynomial P){
	real_polynomial DP;
	DP=D(P);	// derivative
	double t,u;
	t=0.0;
	u=1.0;
	while(EVAL(P,t)*EVAL(P,u)>= 0.0){	// double u until the root is between t and u
		u=u*2.0;
	};
	while(EVAL(P,u)>0.0001){	// don't need much accuracy here.
//		cout << "u " << u << " P(u) " << EVAL(P,u) << "\n";
		u=u-EVAL(P,u)/EVAL(DP,u);
	};
	return(u);
};

real_polynomial J_T_real_poly(polynomial P){
	int i;
	real_polynomial Q;
	Q.a.resize(0);
	for(i=0;i<=P.degree();i++){
		Q.a.push_back(abs(P.a[i]));
	};
	Q.a[0]=-Q.a[0];
	return(Q);
};

complex<double> J_T_find_root(polynomial P){
/* computes roots by Jenkins-Traub algorithm. See
	 Jenkins, M. A. and Traub, J. F. (1970), A Three-Stage Variables-Shift Iteration 
	 for Polynomial Zeros and Its Relation to Generalized Rayleigh Iteration, 
	 Numer. Math. 14, 252–263.	*/

	polynomial Q,H,J;
	real_polynomial K;
	complex<double> s,z,I,t,tt;
	int i;

	I.real()=0.0;
	I.imag()=1.0;
	
	Q=P*(1.0/P.a[P.degree()]);	// make monic
	
	H=Q.D()*(1.0/(double) P.degree());	// initialize H_0 = P'
	s=0.0;
	
	// stage 1
	for(i=0;i<5;i++){	// M=5?
		J=Q-H*(Q.EVAL(0)/H.EVAL(0));
		H=J/monomial(1.0,1);			
	};

	// stage 2
	K=J_T_real_poly(Q);
	s=find_real_positive_simple_root(K);	//  Newton's method. But it's OK because K is real.
	
	s=s*exp(I*(double) (100.0*rand() / RAND_MAX));	// give it a random argument
	t=s;
	for(i=0;i<9;i++){	// L=14?
		J=Q-H*(Q.EVAL(s)/H.EVAL(s));
		H=J/(monomial(1.0,1)-monomial(s,0));
	};	
	// stage 3
	i=14;
	while(norm(Q.EVAL(s))>0.00000000000001){
		J=Q-H*(Q.EVAL(s)/H.EVAL(s));
		H=J/(monomial(1.0,1)-monomial(s,0));
//		H=H*(1.0/H.a[H.degree()]);
		z=H.a[H.degree()];
		s=s-(z*Q.EVAL(s)/H.EVAL(s));
		i++;
		if(i>100){
			cout << "numerical error; more than 100 steps to find root with J-T.\n";
			return(s);	// this is the best we got
		};
	};

	return(s);
};

vector<complex<double> > J_T_find_roots(polynomial P){
	polynomial Q;
	vector<complex<double> > r;
	complex<double> z;
	Q=P;
	r.clear();
	while(Q.degree()>0){
		z=J_T_find_root(Q);
		r.push_back(z);
		Q=Q/(monomial(1.0,1)-monomial(z,0));
	};
	return(r);
};
