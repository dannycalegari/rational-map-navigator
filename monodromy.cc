/* monodromy.cc		algorithm to determine braid sequence to adjust monodromy to desired value */

struct transposition{	// interchange i with j (order doesn't matter)
	int i,j;	 // ASSUME i<j for concreteness!!
};

void write_transposition_sequence(vector<transposition> T){
	int i;
	for(i=0;i<(int) T.size();i++){
		cout << "(" << T[i].i << "," << T[i].j << ") ";
	};
	cout << "\n";
};

void enforce_order(transposition &t){	// switches i and j if necessary so i<j
	int k;
	if(t.i>t.j){
		k=t.i;
		t.i=t.j;
		t.j=k;
	};
};

class permutation{
	public:
		vector<int> i;	
		int operator()(int j){
			return(i[j]);
		};
};

void write_permutation(permutation P){
	int i;
	for(i=0;i<(int) P.i.size();i++){
		cout << i << "->" << P(i) << " ";
	};
	cout << "\n";
};

struct braid{	// move critical value i j units to the right (j can be negative)
	// braid under if k=-1, and over if k=1.
	int i,j;
	bool over;
};

void write_braid(braid b){
	cout << "start at " << b.i << ", move " << b.j << ", over=" << b.over << "\n";
};

void write_braid_sequence(vector<braid> B){
	int i;
	for(i=0;i<(int) B.size();i++){
		cout << "braid " << i << " ";
		write_braid(B[i]);
	};
};

void prune_braid_sequence(vector<braid> &B){	// eliminate braids which do nothing (i.e. j=0)
	int i;
	for(i=0;i<(int) B.size();i++){
		if(B[i].j==0){
			B.erase(B.begin()+i);
			i--;
		};
	};
};

void conjugate(vector<transposition> &T, permutation P){	// changes T to T^P
	int i;
//	cout << "performing conjugation by permutation \n";
	write_permutation(P);
	for(i=0;i<(int) T.size();i++){
		T[i].i=P(T[i].i);
		T[i].j=P(T[i].j);
		enforce_order(T[i]);
	};
	write_transposition_sequence(T);
};

transposition conjugate(transposition t, transposition s){	// returns t^s
	transposition r;
	r=t;
	if(t.i==s.i){
		r.i=s.j;
	};
	if(t.i==s.j){
		r.i=s.i;
	};
	if(t.j==s.i){
		r.j=s.j;
	};
	if(t.j==s.j){
		r.j=s.i;
	};	
	enforce_order(r);
	return(r);
};

void operate(vector<transposition> &T, braid &b){	// change a transposition sequence by performing a braid
	int n;
	transposition s;
//	cout << "operating on S with braid \n";
	write_braid(b);
	while(b.j!=0){
		if(b.j>0){	// move to the right
			n=(b.i+1) % (int) T.size();
		} else {
			n=((int) T.size() + b.i-1) % (int) T.size();
		};
		if(b.over){
			s=T[b.i];
			T[b.i]=T[n];
			T[n]=conjugate(s,T[n]);
		} else {
			s=T[b.i];
			T[b.i]=conjugate(T[n],T[b.i]);
			T[n]=s;
		};
		b.i=n;
		if(b.j>0){
			b.j--;
		} else {
			b.j++;
		};
	};
	write_transposition_sequence(T);
};



int find_matching_transposition(vector<transposition> T, int start_from, int i, int LEGi, int j, int LEGj){
	/* T is a transposition sequence of degree d.
	This function finds the smallest k>=i such that T[k] matches the pattern (LEGi i), (LEGj j)
	where LEGi and LEGj are one of {-1,0,1}. Here -1 stands for <, 0 stands for =, 1 stands for >
	*/
	int k;
	int n,m;
	
	for(k=start_from;k<(int) T.size();k++){
		n=T[k].i;
		m=T[k].j;
	//	cout << "examining location " << k << " j=" << n << " k=" << m << "\n";
	//	cout << (n>i)-(i>n) << " " << LEGi << " " << (m>j)-(j>m) << " " << LEGj << "\n";
		if((n>i)-(i>n)==LEGi && (m>j)-(j>m)==LEGj){	// cute formula for signum(n-i) and signum(m-j)
			return(k);
			break;
		};
	};
	return(-1);	// indicates we couldn't find matching transposition
};

int find_partial_matching_transposition(vector<transposition> T, int start_from, int  i, int LEGi){
	int k;
	int n,m;
	
	for(k=start_from;k<(int) T.size();k++){
		n=T[k].i;
		m=T[k].j;
		if((n>i)-(i>n)==LEGi || (m>i)-(i>m)==LEGi){	// cute formula for signum(n-i) and signum(m-i)
			return(k);
		};
	};
	return(-1);	// indicates we couldn't find matching transposition
};


vector<braid> compute_reset_sequence(int d, vector<transposition> T){
	/* assuming T is a legal transposition sequence of degree d, returns a vector
	of braids which transform T to the base transposition sequence.
	
	Base transposition sequence is (d-2,d-1) (d-3,d-2) . . . (0,1) (0,1) . . . (d-2,d-1)
	*/
	vector<transposition> S;
	vector<braid> B;
	braid b;
	permutation P;
	vector<int> L;
	
	B.clear();	// B initialized to empty braid
	S=T;	// initial value of S equal to T.
//	cout << "computing reset sequence. current sequence is\n";
	write_transposition_sequence(S);
	
	int i,j,k,l,m,i_j,i_jj,where;
	for(i=0;i<d-1;i++){
//	cout << "adjusting level " << i << "\n";
	// S is in the form (i-1,i) . . . (0,1) (0,1) . . . (i-1,1) *
		if(i==0){	// special case
		//	cout << "i=0; special case \n";
			j=S[0].i;
			k=S[0].j;
			P.i.clear();
			for(l=0;l<d;l++){
				P.i.push_back(l);	// initialize P to identity
			};
			P.i[0]=j;
			P.i[j]=0;
			conjugate(S,P);
			P.i.clear();
			for(l=0;l<d;l++){
				P.i.push_back(l);	// initialize P to identity
			};			
			P.i[1]=S[0].j;
			P.i[S[0].j]=1;
			conjugate(S,P);
		} else {
	//		cout << "at location " << 2*i << " looking for j<=" << i << ", k>" << i << "\n";
			l=find_matching_transposition(S,2*i,i+1,-1,i,1);	
	//		cout << "found it at location " << l << "\n";
	// l is index of first (j,k) in * with j<i+1 and k>i
			if(l>2*i){	
			// move it to the left
			//	cout << "pushing it to the left.\n";
				b.i=l;
				b.j=2*i-l;
				b.over=false;
				B.push_back(b);
				operate(S,b);	// resets b
			};
			j=S[2*i].i;
			k=S[2*i].j;
	// S is in the form (i-1,i) . . . (0,1) (0,1) . . . (i-1,i) (j,k) *
			for(l=0;l<j+1;l++){	// rotate j+1 times
			//	cout << "performing rotation " << l << "\n";
				b.i=i;
				b.j=i-1;
				b.over=true;
				B.push_back(b);
				operate(S,b);
				b.i=i-1;
				b.j=1-i;
				b.over=true;
				B.push_back(b);
				operate(S,b);
			};
	// S is in the form (j-1,j) . . . (j,j+1) (j,j+1) . . . (j-1,j) (j,k) *
			P.i.clear();
			for(l=0;l<d;l++){
				P.i.push_back(l);	// initialize P to identity
			};
			for(l=0;l<=i;l++){
				P.i[l]=((l+i-j) % (i+1));	// P takes l to l+i-j % i for 0<=l<=i
			};
			P.i[i+1]=k;	// P takes i+1 to k
			P.i[k]=i+1;	// P takes k to i+1
			conjugate(S,P);
		};
		
	// S is in the form (i-1,i) . . . (0,1) (0,1) . . . (i-1,i) (i,i+1) *
		L.resize(0);
	// look for a sequence of the form (i+1,i_0) * (i_0,i_1) * . . * (i_j,i)
		i_j=i+1;
		where=2*i;
		while(i_j!=i){
		//	cout << "at location " << where << " looking for " << i_j << "\n";
			l=find_partial_matching_transposition(S,where+1,i_j,0);
		//	cout << "found match at location " << l << "\n";
			// S(l) = (i_j,i_{j+1}) or (i_{j+1},i_j)
			if(S[l].i==i_j){
				i_jj=S[l].j;
				where=l;
			} else if (S[l].j==i_j){
				i_jj=S[l].i;
				where=l;
			} else {
			//	cout << "ERROR; couldn't find partial matching transposition.\n";
				break;
			};
			if(i_jj==i+1){
				L.resize(0);	// initialize vector to 0 (want rightmost sequence)
			} else {
				L.push_back(l);
				i_j=i_jj;
				where=l;
			};
		};
	// sequence (i+1,i_0) * (i_0,i_1) * . . * (i_j,i) occurs at L[0], L[1], . . L[L.size()-1]
	// take L[0] and slide it to the right, under each * and over each L[k]
	//	cout << "L sequence is ";
	//	for(l=0;l<(int) L.size();l++){
	//		cout << L[l] << " ";
	//	};
	//	cout << "\n";
		where=L[0];
		for(l=1;l<(int) L.size();l++){
			m=L[l];
			b.i=where;
			b.j=m-where-1;
			b.over=false;
			B.push_back(b);
			operate(S,b);
			where=m-1;
			b.i=where;
			b.j=1;
			b.over=true;
			B.push_back(b);
			operate(S,b);
			where++;
		};
	//	cout << "where = " << where << "\n";
		b.i=where;
		b.j=2*d-3-where;
		b.over=false;
		B.push_back(b);
		operate(S,b);
	// S is in the form (i-1,i) . . . (0,1) (0,1) . . . (i-1,i) (i,i+1) * (i+1,i)
		b.i=2*d-3;
		b.j=3-2*d;
		b.over=true;
		B.push_back(b);
		operate(S,b);
	// S is in the form (i,i+1) (i-1,i) . . . (0,1) (0,1) . . . (i-1,i) (i,i+1) * 
	};
	prune_braid_sequence(B);
	cout << "reset sequence is \n";
	write_braid_sequence(B);
	return(B);
};