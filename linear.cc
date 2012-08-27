/* linear.cc linear algebra operations */

int closest_entry(vector<complex<double> > V, complex<double> z){	
	// returns i minimizing abs(V[i]-z)
	int i,j;
	double t;
	i=0;
	t=abs(V[0]-z);
	for(j=0;j<(int) V.size();j++){
		if(abs(V[j]-z)<t){
			i=j;
			t=abs(V[j]-z);
		};
	};
	return(i);
};
	
vector<complex<double> >  mul_matrix(vector< vector<complex<double> > > N, vector<complex<double> > v){
	// returns vector Nv
	vector<complex<double> > w;
	complex<double> z;
	int i,j;
	w.resize(0);
	for(i=0;i<(int) N[0].size();i++){
		z=0.0;
		for(j=0;j<(int) N.size();j++){
			z=z+N[j][i]*v[j];
		};
		w.push_back(z);
	};
	return(w);
};

vector<complex<double> >  invert_matrix(vector< vector<complex<double> > > N, vector<complex<double> > v){
	/* 	Given an n x (n+r) matrix N of full rank n, and a column vector v of length n with a single nonzero entry,
	we want to find something in N^{-1}(v). We do this by column operations. 
	
	Basically, we start with a copy M of N and U initialized to the (n+r) x (n+r) identity matrix. Then
	we perform an identical sequence of column operations to M and to U so that at the end,
	the first column of M is equal to v. Since NU=M stays constant throughout, the first 
	column of U is our desired answer.	
	
	Note that this implementation ASSUMES, but does NOT check, that
	N is square of full rank, and has the same size as v. */	

	// returns vector U[0] such that NU[0] = v. I.e. computes u = N^{-1}v.
	// N[j] is column j, N[j][i] is column j row i.
	
	int rows, cols;
	int j,i,k,m;
	double t;
	vector<complex<double> > w;	// swap column
	complex<double> z;	// swap complex
	complex<double> l;	
	
	vector< vector<complex<double> > > M;
	M=N;	// initialize M to equal N
	
	rows=M[0].size();
	cols=M.size();		// note! cols>rows!
	
	vector< vector<complex<double> > > U;	// matrix
	for(j=0;j<cols;j++){	// initialize U to the identity matrix
		w.resize(0);
		for(i=0;i<cols;i++){
			w.push_back(0.0);
		};
		w[j]=1.0;
		U.push_back(w);
	};
	w.resize(0);
	
	// RULE: as we adjust M, we adjust U so that always  NU = M

	/* Step 1: make M lower-diagonal. For each row i in turn, we permute columns
	j >= i so that the i,i entry is at least as big as the i,j entry (this reduces round off
	errors, since we never need to divide by a very small number; see
	
	J. von Neumann, H.H. Goldstine, Numerical Inverting of Matrices of High Order, 
	Bull. Amer. Math. Soc., Vol. 53, No. 11, pp. 1021-1099, 1947.
	
	Then we subtract off multiples of column i from columns j so that i,j entry is zero for j>i. */
	
	for(i=0;i<rows;i++){	// for each row i in turn
	
		// find j>=i maximizing |M[j][i]|
		k=i;
		t=abs(M[k][i]);
		for(j=i;j<cols;j++){
			if(abs(M[j][i])>t){
				k=j;
				t=abs(M[j][i]);
			};
		};
	// swap columns k and i
		w.resize(0);
		for(m=0;m<rows;m++){
			w.push_back(M[i][m]);
			M[i][m]=M[k][m];
			M[k][m]=w[m];
		};
	// keep U in sync
		w.resize(0);
		for(m=0;m<cols;m++){
			w.push_back(U[i][m]);
			U[i][m]=U[k][m];
			U[k][m]=w[m];
		};
	
		for(j=i+1;j<cols;j++){	// add multiples of M[i] to M[j] for j>i so M[j][i]=0
			l=M[j][i]/M[i][i];
			// subtract l.M[i] from M[j]
			for(k=0;k<rows;k++){
				M[j][k]=M[j][k]-l*M[i][k];
			};
			// keep u in sync
			for(k=0;k<cols;k++){
				U[j][k]=U[j][k]-l*U[i][k];
			};
		};
	};

	/* At this point M is lower-diagonal. 
	Step 2: make first column of M equal to v. For each row i we add a multiple of column i
	to column 0 so that the i entry of column 0 equals to i entry of v. */
	
	for(i=0;i<rows;i++){	// for each row in turn
		l=(v[i]-M[0][i])/M[i][i];		// find l such that l.M[i][i] + M[0][i] = v[i]
	// add l times column i to column 0
		for(k=0;k<rows;k++){
			M[0][k]=M[0][k]+l*M[i][k];
		};
	// keep u in sync
		for(k=0;k<cols;k++){
			U[0][k]=U[0][k]+l*U[i][k];
		};
	};
	// at this point first vector of M should be v. So first vector of U is what we want.
	
	return(U[0]);
};
