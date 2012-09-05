/* roots.cc	functions for polynomials and roots */

// make monic polynomial from roots

poly vtop(cvec R, cpx m){	// makes polynomial P with roots R and leading coefficient m
	poly P;
	P=ctop(m);	// P=m
	int i;
	for(i=0;i<(int) R.size();i++){
		P=P*factor(R[i]);	// P=P*(z-R[i])
	};
	return(P);
};

// real polynomial functions for brief use in Jenkins-Traub algorithm

int deg(rpoly P){	// degree of real poly
	return((int) P.size()-1);
};

rpoly D(rpoly P){	// derivative of real poly
	rpoly R;
	int i;
	R.resize(0);
	if(deg(P)==0){
		R.push_back(0);
	} else {
		for(i=1;i<=deg(P);i++){
			R.push_back(P[i] * (double) i);		
		};
	};
	return(R);
};

double EVAL(rpoly P, double t){	// P(w) 
	int i;  
	double u;
	u=0.0;
	for(i=(int) P.size()-1;i>=0;i--){
		u=(u*t) + P[i];
	};
	return(u);
};

double positive_root(rpoly P){
	rpoly DP;
	DP=D(P);	// derivative
	double t,u;
	t=0.0;
	u=1.0;
	while(EVAL(P,t)*EVAL(P,u)>= 0.0){	// double u until the root is between t and u
		u=u*2.0;
	};
	while(EVAL(P,u)>0.0001){	// don't need much accuracy here.
		u=u-EVAL(P,u)/EVAL(DP,u);
	};
	return(u);
};

rpoly J_T_real_poly(poly P){	// special real polynomial for J-T algorithm
	int i;
	rpoly Q;
	Q.resize(0);
	for(i=0;i<=deg(P);i++){
		Q.push_back(abs(P[i]));
	};
	Q[0]=-Q[0];
	return(Q);
};

// Jenkins-Traub algorithm to find roots of polynomial

cpx find_root(poly P, cpx seed, double accuracy){

	/* computes roots by Jenkins-Traub algorithm. See
	 Jenkins, M. A. and Traub, J. F. (1970), A Three-Stage Variables-Shift Iteration 
	 for Polynomial Zeros and Its Relation to Generalized Rayleigh Iteration, 
	 Numer. Math. 14, 252–263.	*/
	poly Q,H,J;
	rpoly K;
	cpx s,ss,z;
	int i;
	cpx t,tt,ttt;
	
	// first try Newton's method with the seed
	s=seed;
	for(i=0;i<10;i++){
		s=s-EVAL(P,s)/DEVAL(P,s);
		if(norm(EVAL(P,s))<accuracy){
			return(s);
		};
	};
	
	Q=P*(1.0/P[deg(P)]);	// make monic
	H=D(Q)*(1.0/(double) deg(P));	// initialize H_0 = P'/deg(P)
	s=0.0;
	// stage 1
	for(i=0;i<5;i++){	// M=5?
		z=EVAL(Q,0.0)/EVAL(H,0.0);
		J=Q-(H*z);
		H=J/(factor(0.0));		
	};

	// stage 2
	K=J_T_real_poly(Q);
	ss=positive_root(K);	//  Newton's method. But it's OK because K is real.

	ss=ss*exp(I*0.85521133347722);	// 49 degrees
	t=0.0;
	tt=0.0;
	ttt=0.0;
	while(1){
		ss=ss*exp(I*1.64060949687467);
		s=ss;
		for(i=0;i<5;i++){
			ttt=tt;
			tt=t;
			J=Q-H*(EVAL(Q,s)/EVAL(H,s));
			H=J/(monomial(1.0,1)-monomial(s,0));	
			t=s-(H[deg(H)]*EVAL(Q,s)/EVAL(H,s));
		};
		if(norm(tt-ttt)<=norm(ttt)/4.0 && norm(t-tt)<=norm(tt)/4.0){
			i=10;
	// stage 3
			s=s-(H[deg(H)]*EVAL(Q,s)/EVAL(H,s));
			while(norm(EVAL(Q,s))>accuracy){
				J=Q-H*(EVAL(Q,s)/EVAL(H,s));
				H=J/factor(s);
				H=H*(1.0/H[deg(H)]);
				z=H[deg(H)];
				s=s-(z*EVAL(Q,s)/EVAL(H,s));
				i++;
				if(i%40==0){
					cout << "numerical error; more than " << i << " steps to find root with J-T.\n";
					cout.flush();
					cout << "Q is ";
					write(Q);
					cout << "P is ";
					write(P);
					cout << "seed is " << seed << "\n";
					cout << "norm(Q(s)) is " << norm(EVAL(Q,s)) << "\n";
					assert(1==0);
				};
			};
			return(s);
		};
	};
};

cvec factorize(poly P, cvec seed, double accuracy){	// find roots of polynomial P; can take seed of vectors as initial guess
	cvec roots;
	int i;
	cpx v;
	poly Q;
	roots.resize(0);
	Q=P;
	for(i=0;i<deg(P);i++){
		if((int) seed.size()>i){
			v=find_root(Q,seed[i],accuracy);
		} else {
			v=find_root(Q,0.0,accuracy);
		};
		roots.push_back(v);
		Q=Q/factor(v);
	};
	return(roots);
};