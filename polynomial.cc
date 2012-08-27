/*	functions.cc	functions of complex parameters */

class polynomial{
	public:
		vector<complex<double> > a;		// coefficients
		vector<complex<double> > r;		// roots
		complex<double> m;				// multiplier

		complex<double> operator()(complex<double> z){		// want to be able to write P(z)
			int i;	
			complex<double> w (0.0,0.0);
			for(i=(int) a.size()-1;i>=0;i--){
				w=(w*z) + a[i];
			};
			return(w);
		};	
		complex<double> EVAL(complex<double> z){	// want to be able to write EVAL(z) in internal functions
			int i;	
			complex<double> w (0.0,0.0);
			for(i=(int) a.size()-1;i>=0;i--){
				w=(w*z) + a[i];
			};
			return(w);		
		};
		int degree();					// degree of polynomial
		void compute_coefficients();	// determine coefficients a[*] from roots r[*] and multiplier m
		polynomial D();					// derivative
		complex<double> find_root();	// find a root by Newton's method
		complex<double> find_nearby_root(complex<double>);	// find a root by Newton's method, with initial seed value; 
			// note: "find_nearby_root" is surprisingly useless for tracking critical points/values
		complex<double> closest_root(complex<double>);		// compares roots to find closest one to seed

		void compute_roots();			// determine roots r[*] from coefficients a[*]
		void write();
		
		~polynomial(){ };	
};

int polynomial::degree(){
	return(a.size()-1);
};

polynomial monomial(complex<double> c, int d){	// build a polynomial c.z^d
	polynomial P;
	int i;
	for(i=0;i<d;i++){
		P.a.push_back(0.0);
	};
	P.a.push_back(c);
	return(P);
};

void polynomial::write(){
	int i;
	for(i=0;i<=degree();i++){
		cout << a[i] << ".z^" << i;
		if(i<degree()){
			cout << " + ";
		} else {
	//		cout << "\n";
		};
	};
};

complex<double> operator^(complex<double> z, int n){	// returns z^n
	if(n==0){
		return(1.0);
	} else {
		return(z*(z^(n-1)));
	};
};

polynomial operator+(polynomial P, polynomial Q){
	polynomial R;
	R=P;
	int i;
	for(i=0;i<=Q.degree();i++){
		if(R.degree()<Q.degree()){
			R.a.push_back(Q.a[i]);
		} else {
			R.a[i]=R.a[i]+Q.a[i];
		};
	};
	return(R);
};

polynomial operator-(polynomial P, polynomial Q){
	polynomial R;
	R=P;
	int i;
	for(i=0;i<=Q.degree();i++){
		if(R.degree()<Q.degree()){
			R.a.push_back(-Q.a[i]);
		} else {
			R.a[i]=R.a[i]-Q.a[i];
		};
	};
	return(R);
};

polynomial operator*(polynomial P, polynomial Q){
	polynomial R;
	R.a.resize(0);
	complex<double> z;
	int i,j;
	for(i=0;i<=P.degree()+Q.degree();i++){
		z=0.0;
		for(j=0;j<=P.degree();j++){
			if(i-j>=0 && i-j<=Q.degree()){
				z=z+(P.a[j]*Q.a[i-j]);
			};
		};
		R.a.push_back(z);
	};
	return(R);
};

polynomial operator/(polynomial P, polynomial Q){	// returns R so degree(P-R*Q) < degree(P)
	polynomial R,PP;
	if(P.degree()<Q.degree()){
		R.a.push_back(0.0);
		return(R);
	} else {
		R=monomial(P.a[P.degree()]/Q.a[Q.degree()],P.degree()-Q.degree());
		PP=P-R*Q;
		PP.a.pop_back();
		return(R+PP/Q);
	};
};

polynomial operator%(polynomial P, polynomial Q){ // returns R so P=SQ+R with degree(R) < degree(Q)
	polynomial R;
	int i;
	R=P-(Q*(P/Q));
	for(i=0;i<=P.degree()-Q.degree();i++){
		R.a.pop_back();
	};
	return(R);
};

polynomial operator*(polynomial P, complex<double> z){
	polynomial Q;
	Q=monomial(z,0);
	return(P*Q);
};

void polynomial::compute_coefficients(){		// determine coefficients a[*] from roots r[*] and multiplier m
	polynomial P,L;
	int i;
	P.a.resize(0);	// first, clear coefficients
	P.a.push_back(1.0);	// initialize to the constant polynomial 1
	for(i=0;i<(int) r.size();i++){	// for each root
		L.a.resize(0);
		L.a.push_back(-r[i]);
		L.a.push_back(1.0);		// L = (z-r[i])
		P=P*L;
	};
	a.resize(0);
	for(i=0;i<(int) P.a.size();i++){
		a.push_back(m*P.a[i]);
	};
};

polynomial make_polynomial(vector<complex<double> > r, complex<double> m){	// returns a polynomial with roots r and multiplier m
	polynomial P;
	int i;
	P.m=m;
//	P.r.resize(0);
	for(i=0;i<(int) r.size();i++){
		P.r.push_back(r[i]);
	};
	P.compute_coefficients();
	return(P);
};

polynomial polynomial::D(){		// derivative of polynomial
	polynomial Q;
	int i;
	Q.a.resize(0);
	for(i=1;i<=degree();i++){
		Q.a.push_back((complex<double>) i*a[i]);
	};
	return(Q);
};

polynomial D(polynomial P){		// derivative of polynomial
	polynomial Q;
	int i;
	Q.a.resize(0);
	for(i=1;i<=P.degree();i++){
		Q.a.push_back((complex<double>) i*P.a[i]);
	};
	return(Q);
};
	
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
			z=z-(EVAL(z))/(D().EVAL(z));	// standard adjustment
			i++;
			if(i>100){		// if it hasn't found a root quickly, pick new initial value
				i=0;
				real(z)=(double) (100.0*rand() / RAND_MAX)-50.0;
				imag(z)=(double) (100.0*rand() / RAND_MAX)-50.0;
			};
		};
		return(z);
	};
};

complex<double> polynomial::find_nearby_root(complex<double> seed){	// finds a root by Newton's method
	complex<double> z;
	int i;
//	cout << "initial value " << seed.real() << " + " << seed.imag() << "i\n";
//	cout << "EVAL is " << EVAL(seed).real() << " + " << EVAL(seed).imag() << "i\n";
	if(degree()==0){
		cout << "degree 0 has no roots! \n";
		return(0.0);
	} else {			// find a root by Newton's method
		z=seed;
		while(norm(EVAL(z))>0.0000000000000000001){
			z=z-(EVAL(z))/(D().EVAL(z));	// standard adjustment
			i++;
			if(i>100){		// if it hasn't found a root quickly, pick new initial value
				cout << "\n";
				cout << "	ERROR: couldn't find a nearby root! picking random seed!\n";
				cout << "\n";
				i=0;
				real(z)=(double) (100.0*rand() / RAND_MAX)-50.0;
				imag(z)=(double) (100.0*rand() / RAND_MAX)-50.0;
			};
		};
		return(z);
	};
};

complex<double> polynomial::closest_root(complex<double> z){	// compares roots to find closest one to seed
	complex<double> w;
	int i, closest;
	double closest_dist;
	closest=0;
	closest_dist=abs(r[0]-z);
	for(i=0;i<(int) r.size();i++){
		if(abs(r[i]-z)<closest_dist){
			closest=i;
			closest_dist=abs(r[closest]-z);
		};
	};
	return(r[closest]);
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


polynomial Wronskian(polynomial P, polynomial Q){	// P'Q - Q'P
	polynomial W;
	W=D(P)*Q - D(Q)*P;
	W.a.pop_back();		// biggest coefficient will be zero
	return(W);
};

complex<double> stereo_point(complex<double> z){
	complex<double> w;
	double d;
	d=abs(z);
	w.real()=z.real()*3.0/(1.0+d);
	w.imag()=z.imag()*3.0/(1.0+d);
	return(w);
};

complex<double> inverse_stereo(complex<double> z){
	complex<double> w;
	double d;
	d=abs(z);
	w.real()=z.real()/(3.0-d);
	w.imag()=z.imag()/(3.0-d);
	return(w);
};
