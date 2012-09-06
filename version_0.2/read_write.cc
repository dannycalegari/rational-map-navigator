/* read_write.cc	*/

void rational_map::read_from_file(ifstream &input_file){		// read data from file
	int d,i;
	cpx u;
	input_file >> d;	// degree
	ZERO.resize(0);
	ZERO.push_back(0.0);	// first zero is at 0

	for(i=0;i<d-1;i++){
		input_file >> u.real() >> u.imag();
		ZERO.push_back(u);	
	};
	
	POLE.resize(0);
	for(i=0;i<d-1;i++){
		input_file >> u.real() >> u.imag();
		POLE.push_back(u);
	};
	// last pole is at infinity
	
	input_file >> M.real() >> M.imag();
	CRIT.resize(0);
	VALS.resize(0);
	P.resize(0);
	Q.resize(0);
	G.user_control=false;
	G.labels_on=true;
	G.doing='U';
	
	compute_P_and_Q();			
	compute_critical_points(0.0000001);		// default accuracy
	compute_critical_values();	
};


void rational_map::write_to_file(ofstream &output_file){		// write data to file
	int d,i;
	d=(int) ZERO.size();
	output_file << d;	// degree
	output_file << "\n\n";
	for(i=1;i<d;i++){
		output_file << ZERO[i].real() << "   " << ZERO[i].imag() << "\n";
	};
	output_file << "\n";
	for(i=0;i<d-1;i++){
		output_file << POLE[i].real() << "   " << POLE[i].imag() << "\n";
	};
	output_file << "\n";

	output_file << M.real() << "   " << M.imag() << "\n";	
};