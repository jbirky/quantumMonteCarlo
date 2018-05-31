// ===================================
// JESSICA BIRKY (A13002163)
// ===================================

#include "setup.cpp"
#include "utils.cpp"
using namespace std;


// ===================================
// DECLARE FUNCTIONS
// ===================================

vector<double> gaussianSample2D();								// draw sample from 2D gaussian 

double energyClassical(vector<double> vec);
double energyQuantum(double n);
vector<double> returnClassicalDist(double T);
vector<double> returnQuantumDist(double T);
vector<double> estimateError(vector<double> vec);
vector<double> removeAutocorrelation(vector<double> vec);
vector<double> expected(string dist_type);

// From utils.cpp
double get_random(double min, double max);						// get random uniform number
double printMatrixR(vector<double> k);
void   saveFileR(vector<double> vec, string save_name);

vector<double> vec_sum(vector<double> m1, vector<double> m2);
vector<double> vec_diff(vector<double> m1, vector<double> m2);
double vec_avg(vector<double> vec);
double vec_std(vector<double> vec);
double dot(vector<double> m1, vector<double> m2); 				// dot product of two real vectors


// ===================================
// RUN SIMULATION
// ===================================

int main() {

	clock_t begin = clock();
	srand(time(0));				// random number seed
	
	// ========================================= 
	// Question 1
	// vector<double> energies = returnClassicalDist(Tm);
	// vector<double> sample = removeAutocorrelation(energies);
	// string e_save = "classical/energies.dat";
	// saveFileR(sample, e_save);

	// string dist_type = "classical";
	// expected(dist_type);


	// ========================================= 
	// Question 2
	vector<double> energies = returnQuantumDist(Tm);
	vector<double> sample = removeAutocorrelation(energies);
	string e_save = "quantum/energies.dat";
	saveFileR(sample, e_save);

	string dist_type = "quantum";
	expected(dist_type);


	// ========================================= 
	// Question 3



	// ========================================= 
	// Question 4


	// =========================================

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / pow(10,6);
	cout<<"Time elapsed: "<<elapsed_secs<<" sec"<<endl;

	return 0;
}


// ===================================
// MONTE CARLO
// ===================================


double energyClassical(vector<double> vec) {

	double energy = m/2 * pow((vec[3]-vec[2]), 2) + m*pow(w,2)/2 * pow((vec[0]+vec[1])/2, 2);

	return energy;
}


double energyQuantum(double n) {
	
	double energy = (n + .5)*h*w;

	return energy;
}


double energyArbitrary(double n) {

	double energy = (pow(n,2) + .5)*h*w;

	return energy;
}


double energyQuantum2D(double n1, double n2) {

	double energy = (n1 + .5)*h*w + (n2 + .5)*h*w;

	return energy;
}


vector<double> returnClassicalDist(double T) {

	vector<double> energies(Ns,0);
	vector<double> vec = {0,0,0,0};
	vector<double> xvec;
	vector<double> pvec;
	vector<double> prop;
	vector<double> randvec;
	double Ecurr;
	double Eprop;
	double acc_prob;
	double BETA = -1/(kb*T);	

	// counters
	int accept = 0;
	int reject = 0;

	for (int i=0; i<Ns; i++) {

		xvec = gaussianSample2D();
		pvec = gaussianSample2D();
		randvec = {xvec[0], xvec[1], pvec[0], pvec[1]};

		prop  = vec_sum(vec, randvec);
		Ecurr = energyClassical(vec);
		Eprop = energyClassical(prop);

		// Metropolis-Hastings algorithm
		if (Eprop < Ecurr) {
			vec = prop;
			energies[i] = Eprop;
			accept += 1;
		} else {
			acc_prob = exp(BETA*(Eprop-Ecurr));
			double randn = get_random(0,1);
			if (randn < acc_prob) {
				vec = prop;
				energies[i] = Eprop;
				accept += 1;
			} else {
				energies[i] = Ecurr;
				reject += 1;
			}
		}

	}
	double accept_rate = ((double) accept)/ ((double) Ns);

	cout<<accept_rate<<" steps accepted"<<endl;

	return energies;
}


vector<double> returnQuantumDist(double T) {

	vector<double> energies(Ns,0);
	int ncurr = 0;
	int nprop;
	int randint;
	double Ecurr;
	double Eprop;
	double acc_prob;
	double BETA = -1/(kb*T);	

	// counters
	int accept = 0;
	int reject = 0;

	for (int i=0; i<Ns; i++) {

		randint = (rand() %3) - 1;
		nprop = ncurr + randint;
		if (nprop < 0) {
			nprop = 0;
		}
		Ecurr = energyQuantum(ncurr);
		Eprop = energyQuantum(nprop);

		// Metropolis-Hastings algorithm
		if (Eprop < Ecurr) {
			ncurr = nprop;
			energies[i] = Eprop;
			accept += 1;
		} else {
			acc_prob = exp(BETA*(Eprop-Ecurr));
			double randn = get_random(0,1);
			if (randn < acc_prob) {
				ncurr = nprop;
				energies[i] = Eprop;
				accept += 1;
			} else {
				energies[i] = Ecurr;
				reject += 1;
			}
		}

	}
	double accept_rate = ((double) accept)/ ((double) Ns);

	cout<<accept_rate<<" steps accepted"<<endl;

	return energies;
}


vector<double> removeAutocorrelation(vector<double> vec) {

	
	vector<double> sample(Nsamp, 0);

	for (int i=0; i<Nsamp; i++) {
		sample[i] = vec[Nburn + i*Nskip];
	}

	return sample;
}


vector<double> estimateError(vector<double> vec) {

	vector<double> sample = removeAutocorrelation(vec);

	double avg = vec_avg(sample);
	double std = vec_std(sample) /sqrt((double) sample.size());

	vector<double> error = {avg, std};

	return error;
}


vector<double> expected(string dist_type) {

	vector<double> eng_exp(Tm,0);
	vector<double> err_exp(Tm,0);
	vector<double> eng2_exp(Tm,0);
	vector<double> err2_exp(Tm,0);

	vector<double> dist;
	vector<double> error;

	for (int t=0; t<Tm; t++) {
		if (dist_type == "classical") {
			dist  = returnClassicalDist(t);
		} else if (dist_type == "quantum") {
			dist  = returnQuantumDist(t);
		} else {
			cout<<dist_type<<" not a valid distribution"<<endl;
		}
		error = estimateError(dist);

		eng_exp[t] = error[0];
		err_exp[t] = error[1];

		// cout<<t<<": "<<error[0]<<" +/- "<<error[1]<<endl;
	}

	string save_eng = dist_type + "/expected_energy.dat";
	string save_err = dist_type + "/expected_error.dat";
	saveFileR(eng_exp, save_eng);
	saveFileR(err_exp, save_err);
}
