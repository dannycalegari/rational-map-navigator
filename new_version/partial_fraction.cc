/* partial_fraction.cc alternate version of rational map */

class partial_fraction{
	// R(z) = CONS + COEF[0]/(z-POLE[0]) + COEF[1]/(z-POLE[1]) + . . . + COEF[d-1]/(z-POLE[d-1])
	
	public:
		cvec POLE;
		cvec COEF;
		cpx CONS;
		
		cvec CRIT;		// critical points
		cvec VALS;		// critical values
		
		void initialize();
		int wind(cpx, double, int);	// how many critical points in square
		int poles_in(cpx, double);	// how many poles in square
		void compute_critical_points(double);	// need to specify accuracy
		void compute_critical_values();
		
		cpx E(cpx);			// value of R(w)
		cvec E(cvec);		// vector of values of R(W[i])
		cpx DE(cpx);		// value of DR(w)
		cvec DE(cvec);		// vector of values of DR(W[i])
};

void partial_fraction::initialize(){
	POLE.resize(0);
	COEF.resize(0);
	CONS=0.0;
	CRIT.resize(0);
	VALS.resize(0);
};

cpx partial_fraction::E(cpx w){
	int i;
	cpx u;
	u=CONS;
	for(i=0;i<(int) POLE.size();i++){
		u=u + COEF[i]/(w-POLE[i]);
	};
	return(u);
};

cvec partial_fraction::E(cvec W){
	int i;
	cvec V;
	V.resize(0);
	for(i=0;i<(int) W.size();i++){
		V.push_back(partial_fraction::E(W[i]));
	};
	return(V);
};

cpx partial_fraction::DE(cpx w){
	int i;
	cpx u,v;
	u=0.0;
	for(i=0;i<(int) POLE.size();i++){
		v=w-POLE[i];
		u=u-COEF[i]/(v*v);
	};
	return(u);
};

cvec partial_fraction::DE(cvec W){
	int i;
	cvec V;
	V.resize(0);
	for(i=0;i<(int) W.size();i++){
		V.push_back(partial_fraction::DE(W[i]));
	};
	return(V);
};

int partial_fraction::wind(cpx u, double t, int step){	// how many critical points in square
								// with center u, and side lengths t
	int i;
	int n;
	cpx w, v, vv;
	double A;
	n=0;
	w=u-t/2.0-t*I/2.0;
	
	v=DE(w);	// initial value
	A=0.0;
	for(i=0;i<step;i++){
		w=w+t/(double) step;
		vv=DE(w);
		A=A+arg(vv/v);
		v=vv;
	};
	for(i=0;i<step;i++){
		w=w+(t*I)/(double) step;
		vv=DE(w);
		A=A+arg(vv/v);
		v=vv;
	};
	for(i=0;i<step;i++){
		w=w-t/(double) step;
		vv=DE(w);
		A=A+arg(vv/v);
		v=vv;
	};
	for(i=0;i<step;i++){
		w=w-(t*I)/(double) step;
		vv=DE(w);
		A=A+arg(vv/v);
		v=vv;
	};
	return(1.01*A/TWOPI);
};

int partial_fraction::poles_in(cpx u, double t){ // how many poles in square
	// with center u, and side lengths t
	int p,i;
	cpx v;
	p=0;
	for(i=0;i<(int) POLE.size();i++){
		v=POLE[i]-u;
		if(abs(v.real())<=t/2.0 && abs(v.imag())<=t/2.0){
			p=p+2;
		};
	};
	return(p);
};

void partial_fraction::compute_critical_points(double accuracy){	
	// function finds critical points by computing the winding number of arg around
	// the perimeter of boxes in a successively refined mesh
	
	double mesh;
	mesh = 256.0;
	cvec L,K;	// potential places critical points might be
	int i,j,k;
	int p,n;
	int step;
	bool add_point;
	step=25;
	L.resize(0);
	L.push_back(0.0);
	// they are in the squares with center L[i] and sides mesh
	while(mesh>accuracy){
	//	cout << "mesh size " << mesh << "; number of boxes " << (int) L.size() << "\n";
		K.resize(0);
		for(i=0;i<(int) L.size();i++){
			n=wind(L[i],mesh,step);
			p=poles_in(L[i],mesh);
	//		cout << "box with center " << L[i] << " and mesh " << mesh << " ";
	//		cout << "n " << n << " p " << p << "\n";
			if((n+p)>0){	// something here?
	//			cout << "found something.\n";
				K.push_back(L[i]);
			};
		};
		if(mesh>2.0*accuracy){
			L.resize(0);
			for(i=0;i<(int) K.size();i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						L.push_back(K[i]-mesh/4.0 - mesh*I/4.0 + ((double) j * mesh/2.0) + ((double) k * mesh*I/2.0));
					};
				};
			};
		};
		mesh=mesh/1.99;
	};
	CRIT.resize(0);
	for(i=0;i<(int) K.size();i++){
		add_point=true;
		if(norm(DE(K[i]))>0.01){	// don't add poles
			add_point=false;
		};
		for(j=0;j<(int) CRIT.size();j++){
			if(norm(K[i]-CRIT[j])<0.000001){	// don't add the same point twice
				add_point=false;
			};
		};
		if(add_point){
			CRIT.push_back(K[i]);
		};
	};
};

void partial_fraction::compute_critical_values(){
	VALS = partial_fraction::E(CRIT);
};

rational_map pftorm(partial_fraction S){
	rational_map R;
	poly P,L;
	int i,j;
	P.resize(0);
	P.push_back(0.0);
	for(i=0;i<(int) S.POLE.size();i++){
		L.resize(0);
		L.push_back(S.COEF[i]);	// L = COEF[i].(z-POLE[0]).(z-POLE[1]) . . (z-POLE[d-1])/(z-POLE[j])
		for(j=0;j<(int) S.POLE.size();j++){
			if(j!=i){
				L=L*factor(S.POLE[j]);	// L = L.(z-POLE[j])
			};
		};
		P=P+L;
	};
	L.resize(0);
	L.push_back(S.CONS);
	for(i=0;i<(int) S.POLE.size();i++){
		L=L*factor(S.POLE[i]);
	};
	P=P+L;
	R.P=P;
	R.M=P[(int) P.size()-1];	// leading coefficient
	R.Q=vtop(S.POLE,1.0);
	R.POLE=S.POLE;
	R.ZERO.resize(0);
	R.ZERO=factorize(P,R.ZERO);
	R.CRIT=S.CRIT;
	R.VALS=S.VALS;
	return(R);
};