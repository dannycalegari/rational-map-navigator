/* rational_map.cc definitions and functions */

struct graphics_state{	// stuff we need to know for graphics
	bool user_control;
	bool labels_on;
};

class rational_map{
	// R(z) = M.(z-ZERO[0]).(z-ZERO[1]). . . (z-ZERO[d-1])/(z-POLE[0]). . .(z-POLE[d-2])
	//		= M.P(z)/Q(z)
	// where ZERO[d-1]=0
	// and R(1)=1, so M=Q(1)/P(1)

	public:
		cvec ZERO;		// zeros	by convention 
		cvec POLE;		// poles
		cpx M;			// multiplier
		cvec CRIT;		// critical points
		cvec VALS;		// critical values
		poly P;			// numerator
		poly Q;			// denominator
		graphics_state G;	// graphics state
		
		void initialize();
		void compute_P_and_Q();			// would we ever compute P but not Q?
		void compute_critical_points(double);	// parameter for accuracy
		void compute_critical_values();
		
		cpx E(cpx);			// value of R(w)
		cvec E(cvec);		// vector of values of R(W[i])
		cpx DE(cpx);		// value of DR(w)
		cvec DE(cvec);		// vector of values of R(W[i])
		cmat JAC();			// dV/d{M,Z,P}
		
		void flow_VALS_to(cvec, double);	// flow VALS in straight line to specific value
		void draw_state();
};

void rational_map::initialize(){
	ZERO.resize(0);
	POLE.resize(0);
	M=1.0;
	CRIT.resize(0);
	VALS.resize(0);
	P.resize(0);
	Q.resize(0);
	G.user_control=false;
	G.labels_on=true;
};

void rational_map::compute_P_and_Q(){
	P=vtop(ZERO,M);
	Q=vtop(POLE,1.0);
};

void rational_map::compute_critical_points(double accuracy){
	CRIT = factorize(Wronskian(P,Q),CRIT,accuracy);
};

void rational_map::compute_critical_values(){
	cvec TEMP_VALS;
	if((int) VALS.size()==0){	// first time?
		VALS = rational_map::E(CRIT);
	} else {
		TEMP_VALS = rational_map::E(CRIT);
		match_nearby(VALS,TEMP_VALS);
	};
};

cpx rational_map::E(cpx w){
	return(EVAL(P,w)/EVAL(Q,w));
};

cvec rational_map::E(cvec W){
	cvec V;
	int i;
	V.resize(0);
	for(i=0;i<(int) W.size();i++){
		V.push_back(rational_map::E(W[i]));
	};
	return(V);
};

cpx rational_map::DE(cpx w){
//	return(EVAL(Wronskian(P,Q),w)/EVAL(Q*Q,w)); faster evaluation below
	cpx u;
	u=EVAL(Q,w);
	return((DEVAL(P,w)*u - EVAL(P,w)*DEVAL(Q,w))/(u*u));
};

cvec rational_map::DE(cvec W){
	cvec V;
	int i;
	V.resize(0);
	for(i=0;i<(int) W.size();i++){
		V.push_back(rational_map::DE(W[i]));
	};
	return(V);
};

cmat rational_map::JAC(){		// dV/d{Z,P}
	/* This function computes the Jacobian matrix dVALS/d{ZERO,POLE}
	
	For each zero z_i define R(z,l):=R(z-l)/(z-z_i) 
	and define W(z,l):=Wronskian(P(z-l),Q(z-z_i))
	Then R(c+d,z_i+e) = R(c,z_i) + d dR/dz(c,z_i) + e dR/dl(c,z_i) + quadratic terms.
    Also, W(c+d,z_i+e) = W(c,z_i) + d dW/dz(c,z_i) + e dW/dl(c,z_i) + quadratic terms.	
	If we want W(c+d,z_i+e) = 0 + quadratic terms, we must have d = -e (dW/dl)/(dW/dz)[c,z_i].
	
	Hence the derivative dv/dl(z_i), where v(l)=R(c(l)), and c(l) is the zero of W(*,l) near c, is
        dR/dl - (dW/dl)*(dR/dz)/(dW/dz) evaluated at c,z_i.
    Explicitly, this is PQR'/W'(c-z_i)^2 - R/(c-z_i).
        
    For each pole p_i we form R(z,l) = R(z-p_i)/(z-l) and W(z,l) = Wr(P(z-p_i),Q(z-l)). We get
    a corresponding formula R/(c-p_i) - PQR'/W'(c-p_i)^2.   */

    int i,j;
    cpx w,u,v,c;
    cvec ROW;
    cmat JAC;
    JAC.resize(0);
    
	for(i=0;i<(int) CRIT.size();i++){	// for each critical point
		ROW.resize(0);
    	u=DDEVAL(P,CRIT[i])*EVAL(Q,CRIT[i])-EVAL(P,CRIT[i])*DDEVAL(Q,CRIT[i]);	// D(W)(CRIT[i]);
    	v=EVAL(P,CRIT[i])*EVAL(Q,CRIT[i])*DE(CRIT[i])/u;
	   	c=E(CRIT[i]);
		for(j=0;j<(int) ZERO.size()-1;j++){	// last ZERO is fixed at 0
			w=v/((CRIT[i]-ZERO[j])*(CRIT[i]-ZERO[j]));
            w=w-c/(CRIT[i]-ZERO[j]);  
            ROW.push_back(w);		
		};
		for(j=0;j<(int) POLE.size();j++){	// last POLE is infinity (omitted)
   			w=v/((CRIT[i]-POLE[j])*(CRIT[i]-POLE[j]));
            w=w-c/(CRIT[i]-POLE[j]); 
            ROW.push_back(-w);			
		};
		JAC.push_back(ROW);
	};
	
	return(JAC);
};

void rational_map::flow_VALS_to(cvec V, double accuracy){	// flow VALS in straight line to V
	cvec K,L;
	cmat J;
	double SPEED;
	int i,j;
	
	L=V-VALS;	// this is the direction we want to move
	j=0;
	while(norm(L)>accuracy){
		if(j%10==0){
			draw_state();
		};
		J=JAC();	// this is the Jacobian
		K=INV(J,L);		// J*K=L
		SPEED=0.01/sqrt(norm(L));
		if(SPEED>1.0){
			SPEED=1.0;
		};
		for(i=0;i<(int) ZERO.size()-1;i++){
			ZERO[i]=ZERO[i]+K[i]*SPEED;
		};
		for(i=0;i<(int) POLE.size();i++){
			POLE[i]=POLE[i]+K[i-1+(int) ZERO.size()]*SPEED;
		};
		compute_P_and_Q();
		compute_critical_points(accuracy);
		compute_critical_values();

		L=V-VALS;	// this is the direction we want to move
		if(j>10000){
			cout << "numerical error; more than 10000 steps to flow.\n";
			cout << "norm at termination " << norm(L) << "\n";
			assert(1==0);
		};
		j++;
	};
	draw_state();
};

void rational_map::draw_state(){
	int i;
	cpx v;
	point p;
	stringstream T;
	
	erase_field();

	p.x=320;
	p.y=320;
	draw_circle(p,300,0xAAAAAA);
	p.x=p.x+620;
	draw_circle(p,300,0xAAAAAA);
	for(i=0;i<(int) ZERO.size();i++){
		p=cpx_to_point(ZERO[i]);
		draw_concentric_circles(p,2,0x550000);
		if(G.labels_on){
			draw_label(p,i,0x550000);
		};
	};
	for(i=0;i<(int) POLE.size();i++){
		p=cpx_to_point(POLE[i]);
		draw_concentric_circles(p,2,0x000055);
		if(G.labels_on){
			draw_label(p,i,0x000055);
		};
	};
	for(i=0;i<(int) CRIT.size();i++){
		p=cpx_to_point(CRIT[i]);
		draw_concentric_circles(p,2,0x005500);
		if(G.labels_on){
			draw_label(p,i,0x005500);
		};
	};
	for(i=0;i<(int) VALS.size();i++){
		p=cpx_to_point(VALS[i]);
		p.x=p.x+620;
		draw_concentric_circles(p,2,0x550055);
		if(G.labels_on){
			draw_label(p,i,0x550055);
		};
	};
	if(G.labels_on){
		T << "labels are on";
	};
	p.x=30;
	p.y=650;
	draw_text(p,T,0x000000);
	T.str("");

	
	XFlush(display);
};
