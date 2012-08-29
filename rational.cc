/************************************
*                                   *
*	rational.cc version 0.01		*
*	July 30 2012					*
*									*
*   uses circle packing to 			*
*	numerically determine a			*
*	rational map from Hurwitz		*
*	data
*									*
*	Copyright Danny Calegari 2012	*
*	Released under the GPL license	*
*									*
************************************/
	

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <complex>
#include <sstream>

#define debug 0
#define verbose 1
#define PI 3.14159265358979
#define TWOPI 6.28318530717959

using namespace std;
/*
using std::cin;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::complex; */

#include "monodromy.cc"
#include "linear.cc"
#include "points.cc"
#include "polynomial.cc"
#include "rational_map.cc"
#include "graphics.cc"



int main(int argc, char *argv[]){	

	rational_map R;
	bool finished;
	int d;
		
	cout << "Welcome to the rational map explorer!\n";
	cout << "Enter degree of rational map:";
	cin >> d;
	cout << "Setting up initial zeros/poles.\n\n";
	
	R.initialize(d);
	
	setup_graphics();
	setup_font();
	finished=false;
	cout << "Zeros are red, Poles are blue, Critical points are green (all on LHS). \n";
	cout << "Critical values are light blue (on RHS). \n\n";
	
	cout << "Select a zero or pole with mouse button. \n";
	cout << "Adjust its values with arrow keys. \n\n";
	
	cout << "Select a critical value with mouse button. \n";
	cout << "To steer critical values to roots of unity, type [s]. (warning: experimental!) \n";
	cout << "To push critical values out/in, type [a]/[z]. \n";
	cout << "To rotate critical values positively/negatively, type [x]/[c]. \n\n";
	
	cout << "To output data to screen, type [o]. \n";
	cout << "To toggle vector field D/N/S/none, type [v]. \n";
	cout << "To toggle integral curves on/off, type [i] (warning: very slow!) \n";
	cout << "Dial [m] for monodromy. \n\n";

	cout << "To quit, type [q]. \n";
	
	while(finished==false){
		graphics_routine(R,finished);
	};

	
	return(0);
}
