/* polynomial.cc basic functions with complex numbers and polynomials */
/* convention: z always stands for a variable, w for a specific number */

// functions on complex numbers

cpx operator^(cpx w, int j){    // returns w^j
    if(j==0){
		return(1.0);
    } else {
		return(w*(w^(j-1)));
	};
};

cpx cross_ratio(cpx a, cpx b, cpx c, cpx d){	// returns cross-ratio. If b=0,c=1,d=infty this is a
        return( ((a-b)*(c-d))/((c-b)*(a-d)) );
};

// unary functions of polynomials

int deg(poly P){	// degree of polynomial
	return((int) P.size()-1);
};

void write(poly P){		// write polynomial to cout
	int i;
	for(i=0;i<=deg(P);i++){
		cout << P[i] << "z^" << i << " ";
		if(i<deg(P)){
			cout << "+ ";
		};
	};
	cout << "\n";
};

poly D(poly P){		// derivative of polynomial
	poly R;
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

// binary functions of polynomials and complex

cpx EVAL(poly P, cpx w){	// P(w) 
	int i;  
	cpx u (0.0,0.0);
	for(i=(int) P.size()-1;i>=0;i--){
		u=(u*w) + P[i];
	};
	return(u);
};

cpx DEVAL(poly P, cpx w){	// DP(w)
	int i;
	cpx u (0.0,0.0);
	for(i=(int) P.size()-1;i>0;i--){
		u=(u*w) + P[i]*(double) i;
	};
	return(u);
};

cpx DDEVAL(poly P, cpx w){	// DDP(w)
	int i;
	cpx u (0.0,0.0);
	for(i=(int) P.size()-1;i>1;i--){
		u=(u*w) + P[i]*(double) i*(double) (i-1);
	};
	return(u);
};


poly operator*(cpx w, poly P){	// w*P
	poly R;
	int i;
	R=P;
	for(i=0;i<(int) R.size();i++){
		R[i]=w*R[i];
	};
	return(R);
};

poly operator*(poly P, cpx w){	// w*P
	return(w*P);
};

// methods to create simple polynomials

poly ctop(cpx z){	// convert complex to polynomial
	poly P;
	P.resize(0);
	P.push_back(z);
	return(P);
};

poly monomial(cpx w, int j){	// polynomial w.z^j
	poly P;
	int i;
	P.resize(0);
	for(i=0;i<j;i++){
		P.push_back(0.0);
	};
	P.push_back(w);
	return(P);
};

poly factor(cpx w){	// polynomial z-w
	poly P;
	P.resize(0);
	P.push_back(-w);
	P.push_back(1);
	return(P);
};

// binary functions of polynomials

poly operator+(poly P, poly Q){	// P+Q 
	poly R;
	R=P;
	int i;
	for(i=0;i<(int) Q.size();i++){
		if(deg(R)<i){
			R.push_back(Q[i]);
		} else {
			R[i]=R[i]+Q[i];
		};
	};
	return(R);
};

poly operator-(poly P, poly Q){	// P+Q 
	poly R;
	R=P;
	int i;
	for(i=0;i<(int) Q.size();i++){
		if(deg(R)<i){
			R.push_back(-Q[i]);
		} else {
			R[i]=R[i]-Q[i];
		};
	};
	return(R);
};

poly operator*(poly P, poly Q){	// P*Q
	// is it faster if deg(P) is small or deg(Q)?
	poly R;
	R.resize(0);
	cpx w;
	int i,j;
	for(i=0;i<=deg(P)+deg(Q);i++){
		w=0.0;
		for(j=0;j<=deg(P);j++){
			if(i-j>=0 && i-j<=deg(Q)){
				w=w+(P[j]*Q[i-j]);
			};
		};
		R.push_back(w);
	};
	return(R);
};

poly operator/(poly P, poly Q){	// returns R so deg(P-R*Q) < deg(Q)
	poly R,PP;
	cvec S;	// R in opposite order
	int i;
	if(deg(P)<deg(Q)){
		R.resize(0);
		R.push_back(0.0);
		return(R);
	} else {
		S.resize(0);
		PP=P;
		while(deg(PP)>=deg(Q)){
			S.push_back(PP[deg(PP)]/Q[deg(Q)]);
			R=monomial(PP[deg(PP)]/Q[deg(Q)],deg(PP)-deg(Q));
			PP=PP-Q*R;
			PP.pop_back();	// clear top coefficient
		};
		R.resize(0);
		for(i=0;i < (int) S.size();i++){
			R.push_back(S[S.size()-1-i]);	// reverse order of coefficients
		};
		return(R);
	};
};

poly operator%(poly P, poly Q){	// returns R so P=S*Q+R with deg(R)<deg(Q)
	poly R;
	R=P-(Q*(P/Q));
	R.resize((int) Q.size() - 1);
	return(R);
};

poly Wronskian(poly P, poly Q){	// P'Q - Q'P
	poly W;
	W=(D(P)*Q) - (P*D(Q));
	if(deg(P)==deg(Q)){
		W.pop_back();	// top coefficient vanishes
	};
	return(W);
};

