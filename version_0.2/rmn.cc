/*	rmn.cc	

	Rational Map Navigator version 0.02

	September 1 2012

	Copyright Danny Calegari

	released under the terms of the GNU GPL
	
*/

// standard libraries to include

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <complex>
#include <sstream>
#include <limits>
#include <ctime>
#include <assert.h>

using namespace std;

// preprocessor definitions

#define debug 	1							//	debug is true
#define verbose 1							//	verbose is true by default
#define PI 		3.14159265358979
#define TWOPI	6.28318530717959
#define ACC		0.00000000001				// accuracy 1.0e-11
#define cpx		complex<double>				// complex number
#define cvec 	vector<cpx >	// vector of complex numbers
#define cmat	vector< vector<cpx > >		// matrix of complex
/* Instead of defining poly with a struct, define it as a
	vector of complex numbers. We have the convention that the
	degree of the polynomial is the size + 1, and poly[j] is the
	coefficient of z^j. */
#define poly	vector<cpx >	// complex polynomial
#define rpoly	vector<double>				// real polynomial

// global constants

cpx I (0.0,1.0);

#include "braid.cc";
#include "points.cc";
#include "vector.cc";
#include "graphics.cc"
#include "polynomial.cc";
#include "roots.cc";
#include "rational_map.cc";
#include "draw.cc";
#include "insert.cc";
#include "magnify.cc";
#include "read_write.cc";
#include "user_interface.cc";
#include "select.cc";

// #include "benchmark.cc";
// #include "partial_fraction.cc";

int main(int argc, char *argv[]){ 
	rational_map R;
	ifstream input_file;
	
	setup_graphics();

	if(argc>1){
		input_file.open(argv[1]);
		R.read_from_file(input_file);
		input_file.close();
	} else {
		R.initialize();
	};
	while(1){
		R.draw_state();
		R.user_interface();
	};
	
	return(0);
};