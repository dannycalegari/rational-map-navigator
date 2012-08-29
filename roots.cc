/* roots.cc function to compute roots of complex polynomial */


polynomial J_T_real_poly(polynomial P){
	int i;
	polynomial Q;
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
	J=J_T_real_poly(Q);
	s=J.find_root();	//  Newton's method. But it's OK because J is real.
	while(s.real()<=0.0){		// find the unique positive real root.
		J=J/(monomial(1.0,1)-monomial(s,0));
		s=J.find_root();
	};
	
	s=s*exp(I*(double) (100.0*rand() / RAND_MAX));	// give it a random argument
	t=s;
	for(i=0;i<9;i++){	// L=14?
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
