/* benchmark.cc benchmark tests to debug and time numerical routines */

void benchmark_test(){	// benchmark tests
	rational_map R;
	cvec V;
	poly W,P,Q;
	cpx u,v,w;
	cmat J;
	int i,j,d,k;
	clock_t start;
	double accuracy;
	
	for(d=3;d<11;d++){
		cout << "\n";
		cout << "benchmarks for rational map of degree " << d << "\n";
		accuracy=1.0;
		for(k=0;k<1;k++){
	
			R.initialize();
			u=0.0;
			for(i=0;i<d-1;i++){
				u.real()=u.real()+(3.0*rand() / RAND_MAX) - 1.5;
				u.imag()=u.imag()+(3.0*rand() / RAND_MAX) - 1.5;
				R.POLE.push_back(u);
				u.real()=u.real()+(3.0*rand() / RAND_MAX) - 1.5;
				u.imag()=u.imag()+(3.0*rand() / RAND_MAX) - 1.5;
				R.ZERO.push_back(u);
			};
			R.ZERO.push_back(0.0);
			R.M=1.0;
			cout << "\n";
			accuracy=accuracy*0.0001;

			cout << "accuracy " << accuracy << "\n";
			cout << "computing critical values for the first time.\n";
			start=clock();
			R.compute_P_and_Q();
			R.compute_critical_points(accuracy);
			R.compute_critical_values();
			cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";

			cout << "computing adjusted critical values.\n";
			start=clock();
			R.POLE[1]=R.POLE[1]+0.001;
			R.compute_P_and_Q();
			R.compute_critical_points(accuracy);
			R.compute_critical_values();
			cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";
	
			start=clock();
			cout << "computing Jacobian.\n";
			J=R.JAC();
			for(i=0;i<(int) J.size();i++){
				for(j=0;j<(int) J[i].size();j++){
				};
			};
			cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";

			
			cout << "computing flow to roots of unity.\n";
			start=clock();
			V.resize(0);
			for(j=0;j<(int) R.VALS.size();j++){
				u=exp(TWOPI*I*(double) j/(double) R.VALS.size());
				V.push_back(u);
			};
			R.flow_VALS_to(V,accuracy);
			cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";
			R.G.user_control=true;
			R.draw_state();
			usleep(1000000);
			R.G.user_control=false;
		};
	};
};