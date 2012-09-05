/*	rmn.cc	

	Rational Map Navigator version 0.02

	September 5 2012

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
#define TWOPI		6.28318530717959
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

#include "points.cc"
#include "vector.cc";
#include "graphics.cc"
#include "polynomial.cc";
#include "roots.cc";
#include "rational_map.cc";
#include "benchmark.cc";
// #include "partial_fraction.cc";

int main(int argc, char *argv[]){ 
	setup_graphics();
	benchmark_test();
	while(1){
	};
	return(0);
};
