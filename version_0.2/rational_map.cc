/* rational_map.cc definitions and functions */

struct graphics_state{	// stuff we need to know for graphics
	bool user_control;	// is it waiting for user prompt?
	bool labels_on;		// are there labels on the points?
	char doing;			// what are we currently doing?	
		// U = user control, F = flow to roots of unity, I = insert zero/pole, H = help
		// M = magnifying glass, S = select Z/P/V
	double distance;	// while flowing, distance to value
	cpx insert_point;	// center of zero/pole pair
	cpx zero_point;	// zero of zero/pole pair
	point magnify_location;	// location of magnifying glass
	char select_type;	// one of 'Z', 'P', 'V' (zero/pole/vals) or 'N' (none)
	int select_index;
	cpx select_location;
};

class rational_map{
	// R(z) = M.(z-ZERO[0]).(z-ZERO[1]). . . (z-ZERO[d-1])/(z-POLE[0]). . .(z-POLE[d-2])
	//		= M.P(z)/Q(z)
	// where ZERO[0]=0
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
		void normalize();	// move zeros/poles close to unit circle
		
		cpx E(cpx);			// value of R(w)
		cvec E(cvec);		// vector of values of R(W[i])
		cpx DE(cpx);		// value of DR(w)
		cvec DE(cvec);		// vector of values of R(W[i])
		cmat CRIT_JAC();		// dC/d{M,Z,P}
		cmat VALS_JAC();			// dV/d{M,Z,P}
		
		void flow_VALS_to(cvec, double);	// flow VALS in straight line to specific value
		void do_braid(braid);
		void draw_state();
		void read_from_file(ifstream &);		// read data from file
		void write_to_file(ofstream &);		// write data to file
		void user_interface();			// top-level user interaction routine
		void insert_zp();			// insert zero/pole pair
		void magnify();				// magnifying glass
		void select_and_adjust();	// select and adjust location of Z/P/V
};

void rational_map::initialize(){
	ZERO.resize(0);
	ZERO.push_back(0.0);	// first zero is fixed at 0
	POLE.resize(0);
	M=1.0;
	CRIT.resize(0);
	VALS.resize(0);
	P.resize(0);
	Q.resize(0);
	G.user_control=false;
	G.labels_on=true;
	G.doing='U';
};

void rational_map::compute_P_and_Q(){
	P=vtop(ZERO,M);
	Q=vtop(POLE,1.0);
};

void rational_map::compute_critical_points(double accuracy){
	CRIT = factorize(Wronskian(P,Q),CRIT,accuracy);
};

void rational_map::compute_critical_values(){
	VALS = rational_map::E(CRIT);
		/* tracking is done poorly.

	int i,j;
	cvec TEMP_VALS, TEMP_CRIT;
	if((int) VALS.size()==0){	// first time?
		VALS = rational_map::E(CRIT);
	} else {
		VALS = rational_map::E(CRIT);
		TEMP_VALS = rational_map::E(CRIT);
		match_nearby(VALS,TEMP_VALS);	// permutes TEMP_VALS to be near VALS
		TEMP_CRIT=CRIT;
		for(i=0;i<(int) TEMP_CRIT.size();i++){
			j=closest_match(E(TEMP_CRIT[i]),VALS);
			CRIT[j]=TEMP_CRIT[i];
		};
	};
			*/

};

void rational_map::normalize(){		// move zeros/poles close to unit circle
	double t;
	int i;
	t=1.0;
	for(i=1;i<(int) ZERO.size();i++){
		t=t*abs(ZERO[i]);
	};
	for(i=0;i<(int) POLE.size();i++){
		t=t*abs(POLE[i]);
	};
	if(t>1.0){
		t=0.99;
	} else {
		t=1.01;
	};
   	ZERO=ZERO*t;
    POLE=POLE*t;
    M=M/t;
	compute_P_and_Q();
	compute_critical_points(0.00000000000000001);
	compute_critical_values();
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
//	return(EVAL(Wronskian(P,Q),w)/EVAL(Q*Q,w)); // faster evaluation below
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

cmat rational_map::CRIT_JAC(){		// dC/d{Z,P}
	/* This function computes the Jacobian matrix dCRIT/d{ZERO,POLE}
	
	For each zero z_i this is PQR/W'(c-z_i)^2 
        
    For each pole p_i this is -PQR'/W'(c-p_i)^2.   */

    int i,j;
    cpx w,u,v,c;
    cvec ROW;
    cmat JAC;
    JAC.resize(0);
    
	for(i=0;i<(int) CRIT.size();i++){	// for each critical point
		ROW.resize(0);
//		ROW.push_back(0);	// derivative wrt M
    	u=DDEVAL(P,CRIT[i])*EVAL(Q,CRIT[i])-EVAL(P,CRIT[i])*DDEVAL(Q,CRIT[i]);	// D(W)(CRIT[i]);
    	v=EVAL(P,CRIT[i])*EVAL(Q,CRIT[i])/u;
		for(j=1;j<(int) ZERO.size();j++){	// first ZERO is fixed at 0
			w=v/((CRIT[i]-ZERO[j])*(CRIT[i]-ZERO[j]));
            ROW.push_back(-w);		
		};
		for(j=0;j<(int) POLE.size();j++){	// ``-1th'' POLE is infinity (omitted)
   			w=v/((CRIT[i]-POLE[j])*(CRIT[i]-POLE[j]));
            ROW.push_back(w);			
		};
		JAC.push_back(ROW);
	};
	
	return(JAC);
};

cmat rational_map::VALS_JAC(){		// dV/d{Z,P}
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
//		ROW.push_back(CRIT[i]/M);
    	u=DDEVAL(P,CRIT[i])*EVAL(Q,CRIT[i])-EVAL(P,CRIT[i])*DDEVAL(Q,CRIT[i]);	// D(W)(CRIT[i]);
    	v=EVAL(P,CRIT[i])*EVAL(Q,CRIT[i])*DE(CRIT[i])/u;
	   	c=E(CRIT[i]);
		for(j=1;j<(int) ZERO.size();j++){	// first ZERO is fixed at 0
			w=v/((CRIT[i]-ZERO[j])*(CRIT[i]-ZERO[j]));
            w=w-c/(CRIT[i]-ZERO[j]);  
            ROW.push_back(w);		
		};
		for(j=0;j<(int) POLE.size();j++){	// ``-1th'' POLE is infinity (omitted)
   			w=v/((CRIT[i]-POLE[j])*(CRIT[i]-POLE[j]));
            w=w-c/(CRIT[i]-POLE[j]); 
            ROW.push_back(-w);			
		};
		JAC.push_back(ROW);
	};
	
	return(JAC);
};

void rational_map::flow_VALS_to(cvec V, double accuracy){	// flow VALS in straight line to V
	cvec K,LC,L,VV,CC;
	cmat J,JC;
	double SPEED;
	int i,j;
		
	L=V-VALS;	// this is the direction we want to move
	j=0;
	while(norm(L)>accuracy){
		if(j%10==1){
			draw_state();
		};
		J=VALS_JAC();		// dV/d{Z,P}
		JC=CRIT_JAC();	// dC/d{Z,P}
		K=INV(J,L);		// J*K=L
		LC=JC*K;		// how critical points move
		SPEED=0.01/sqrt(norm(L));
		if(SPEED>0.1){
			SPEED=0.1;
		};
		if(SPEED<0.00001){
			SPEED=0.00001;
		};
//		M=M+K[0]*SPEED;
		for(i=1;i<(int) ZERO.size();i++){
			ZERO[i]=ZERO[i]+K[i-1]*SPEED;
		};
		for(i=0;i<(int) POLE.size();i++){
			POLE[i]=POLE[i]+K[i-1+(int) ZERO.size()]*SPEED;
		};
		compute_P_and_Q();
		CRIT=CRIT+(LC*SPEED);	// estimate move first

		VALS=VALS+(L*SPEED);
		compute_critical_points(0.000000000000000001);
		compute_critical_values();
//		normalize();

		L=V-VALS;	// this is the direction we want to move
		G.distance=sqrt(norm(L));
		if(j>10000){
			cout << "numerical error; more than 10000 steps to flow.\n";
			cout << "norm at termination " << norm(L) << "\n";
			cout << "flowing to ";
			vwrite(V);
			assert(1==0);
		};
		if(norm(VALS)>100000){	// something's wrong here
			cout << "norm blow up! dV/d{Z,P} is ";
			J=VALS_JAC();	// this is the Jacobian
			mwrite(J);
			cout << "dC/d{Z,P} is ";
			JC=CRIT_JAC();
			mwrite(JC);
			K=INV(J,L);		// J*K=L
			LC=JC*K;		// how critical points move
			cout << "LC = ";
			vwrite(LC);
			cout << "inverting JAC gives K ";
			vwrite(K);
			cout << "J*K = ";
			vwrite(J*K);
			CC=CRIT;
			VV=VALS;
			SPEED=SPEED/100000.0;
			for(i=1;i<(int) ZERO.size();i++){
				ZERO[i]=ZERO[i]+K[i-1]*SPEED;
			};
			for(i=0;i<(int) POLE.size();i++){
				POLE[i]=POLE[i]+K[i-1+(int) ZERO.size()]*SPEED;
			};
			compute_P_and_Q();
			compute_critical_points(0.000000000000000001);
			compute_critical_values();
			cout << "actual direction of motion is ";
			vwrite((VALS-VV)*(1.0/SPEED));
			cout << "direction L we want to move is ";
			vwrite(L);
			cout << "critical direction of motion is ";
			vwrite((CRIT-CC)*(1.0/SPEED));
			cout << "theoretical direction is ";
			vwrite(LC);
			cout << "sanity check. \n";
			cout << "VALS ";
			vwrite(VALS);
			cout << "E(CRIT) ";
			vwrite(E(CRIT));
			assert(1==0);
		};
		j++;
	};
};

void rational_map::do_braid(braid b){
	cvec V;
	cpx u,v;
	int i,ii,ij,e;
	e=(int) VALS.size();
	u=exp(TWOPI*I/(double) e);	// root of unity
	
	V=VALS;
	if(b.over){		// push V[b.i] in or out
		V[b.i]=V[b.i]*2.0;
	} else {
		V[b.i]=V[b.i]*0.5;
	};
	flow_VALS_to(V,0.001);	// not much accuracy needed
	
	if(b.j>0){		// rotate intervening braids
		for(i=1;i<=b.j;i++){
			ii=(b.i+i)%e;
			V[ii]=V[ii]/u;
		};
	} else {
		for(i=-1;i>=b.j;i--){
			ii=(b.i+i+e)%e;
			V[ii]=V[ii]*u;
		};
	};
	flow_VALS_to(V,0.001);	// not much accuracy needed
	
	if(b.j>0){		// rotate b.i
		for(i=1;i<=b.j;i++){
			V[b.i]=V[b.i]*u;
			flow_VALS_to(V,0.001);	// not much accuracy needed
		};
	} else {
		for(i=-1;i>=b.j;i--){
			V[b.i]=V[b.i]/u;
			flow_VALS_to(V,0.001);	// not much accuracy needed
		};
	};	
	
	V[b.i]=V[b.i]/abs(V[b.i]);		// push V[b.i] out or in
	flow_VALS_to(V,0.001);	// not much accuracy needed


	if(b.j>0){		// reorder indices
		u=VALS[b.i];
		v=CRIT[b.i];
		for(i=1;i<=b.j;i++){
			ii=(b.i+i+e)%e;
			ij=(b.i+i-1+e)%e;
			VALS[ij]=VALS[ii];
			CRIT[ij]=CRIT[ii];
		};
		VALS[(b.i+b.j)%e]=u;
		CRIT[(b.i+b.j)%e]=v;
	} else {
		u=VALS[b.i];
		v=CRIT[b.i];
		for(i=-1;i>=b.j;i--){
			ii=(b.i+i+e)%e;
			ij=(b.i+i+1+e)%e;
			VALS[ij]=VALS[ii];
			CRIT[ij]=CRIT[ii];
		};
		VALS[(b.i+b.j+e)%e]=u;
		CRIT[(b.i+b.j+e)%e]=v;		
	};
	draw_state();
};

