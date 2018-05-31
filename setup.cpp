#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <cblas.h>
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DEFINE PROBLEM CONSTANTS
// ===================================

// Physical constants
double  PI 	= 3.141592652;
double	h	= 1;				// Plancks constant
double	m 	= 1;				// Mass of particle
double	w	= 1; 				// Oscillation frequency
double  kb  = 1;				// Boltzmann constant
int     Tm  = 100;				// Tempertature max of expected plot

// MCMC parameters
int		Ns = pow(10,5);			// Number of steps
int		Nw = 1; 				// Number of walkers
int		Nburn = 1000;			// Burn in 
int 	Nskip = 500; 			// Number steps skipped per iteration
double	DEL_S = 5;				// Width of gaussian step
int 	Nsamp = round((Ns-Nburn)/Nskip); // Number of points kept per sample