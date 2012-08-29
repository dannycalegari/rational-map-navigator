/* roots.cc function to compute roots of complex polynomial */


polynomial J_T_poly(polynomial P){
	int i;
	polynomial Q;
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
	J=J_T_poly(Q);
	s=J.find_root();
	while(s.real()<=0.0){
		J=J/(monomial(1.0,1)-monomial(s,0));
		s=J.find_root();
	};
	s=s*exp(I*(double) (100.0*rand() / RAND_MAX));
	t=s;
	for(i=0;i<9;i++){	// M=5?
		J=Q-H*(Q.EVAL(s)/H.EVAL(s));
		H=J/(monomial(1.0,1)-monomial(s,0));
	};	
	// stage 3
	while(norm(Q.EVAL(s))>0.00000000000001){
		J=Q-H*(Q.EVAL(s)/H.EVAL(s));
		H=J/(monomial(1.0,1)-monomial(s,0));
//		H=H*(1.0/H.a[H.degree()]);
		z=H.a[H.degree()];
		s=s-(z*Q.EVAL(s)/H.EVAL(s));
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

/*
complex<double> polynomial::find_root(){	// finds a root by Newton's method
	complex<double> z;
	int i;
	i=0;
	if(degree()==0){
		cout << "degree 0 has no roots! \n";
		return(0.0);
	} else {			// find a root by Newton's method
		z=0.0;
		while(norm(EVAL(z))>0.0000000000000000001){
			z=z-EVAL(z)/(D().EVAL(z));	// standard adjustment
			i++;
			if(i>100){		// if it hasn't found a root quickly, pick new initial value
				i=0;
				real(z)=(double) (100.0*rand() / RAND_MAX)-50.0;
				imag(z)=(double) (100.0*rand() / RAND_MAX)-50.0;
				cout << "root-finding trouble. \n";
			};
		};
		return(z);
	};
};

void polynomial::compute_roots(){
	polynomial R,S;
	complex<double> z;
	int i;
	r.resize(0);	// clear list of roots
	R.a.resize(0);
	for(i=0;i<(int) a.size();i++){
		R.a.push_back(a[i]);
	};
	while(R.degree()>0){
		z=R.find_root();	// find a root
		r.push_back(z);		// add it to the list of roots
		R=R/(monomial(1.0,1)-monomial(z,0));	// R=R/(z-root)
	};
	m=R.a[0];	// remember multiplier
};

*/