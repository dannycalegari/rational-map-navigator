/* rational_map.cc definitions and functions */

class rational_map{
	public:
		cvec ZERO;		// zeros
		cvec POLE;		// poles
		cpx M;			// multiplier

		cvec CRIT;		// critical points
		cvec VALS;		// critical values

		poly P;			// numerator
		poly Q;			// denominator
		
		void initialize();
		void compute_P_and_Q();	// would we ever compute P but not Q?
		void compute_critical_points();
		void compute_critical_values();
};

void rational_map::initialize(){
	ZERO.resize(0);
	POLE.resize(0);
	M=1.0;
	CRIT.resize(0);
	VALS.resize(0);
	P.resize(0);
	Q.resize(0);
};

void rational_map::compute_P_and_Q()