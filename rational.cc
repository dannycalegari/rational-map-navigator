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
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <complex>

#define debug 0
#define verbose 1
#define PI 3.14159265358979
#define TWOPI 6.28318530717959

using std::cin;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::complex;

#include "polynomial.cc"
#include "rational_map.cc"
#include "graphics.cc"


int main(int argc, char *argv[]){
	vector<complex<double> > roots;
	complex<double> m,z,w,eta;
	polynomial P,Q;
	rational_map R;
	int d,i;
	bool finished;
	
	m=1.0;
	
	cout << "Welcome to the rational map explorer!\n";
	cout << "Enter degree of rational map:";
	cin >> d;
	cout << "Setting up initial zeros/poles.\n";
	
	w.real()=0;
	w.imag()=TWOPI/(2.0*(double) d);
	eta=exp(w);	// 2dth root of unity
	roots.clear();	// initializing roots
	for(i=0;i<d;i++){
		roots.push_back((eta^(2*i)));
		if(i==0){
			roots[0]=roots[0]*1.25;
		};
	};
	R.Zeros=roots;
	roots.clear();
	for(i=0;i<d;i++){
		roots.push_back((eta^(2*i+1)));
		if(i==d/2){
			roots[i]=roots[i]*1.2;
		};
	};
	R.Poles=roots;
	R.M=m;

	R.compute_coefficients();
	R.compute_C_and_V();
	R.ZP='Z';	
	R.ZP_index=0;
	
	setup_graphics();
	finished=false;
	cout << "Zeros are red, Poles are green. \n";
	cout << "Critical points are green, critical values (on RHS) are purple.\n";
	cout << "Select a zero or pole with mouse button. \n";
	cout << "Adjust its values with arrow keys. \n";
	cout << "To output data to screen, type [o]. \n";
	cout << "To quit, type [q]. \n";
	
	while(finished==false){
		graphics_routine(R,finished);
	};
	
	return(0);
}
