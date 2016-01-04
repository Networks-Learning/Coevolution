#include "simulator.cpp"
#include <iostream>
#include <string.h>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
	cout << "Hello Network World!\n";

	double mu_mean, alpha_mean, eta_mean, beta_mean;

	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-N")==0) {
			params.node_count = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-T")==0) {
			params.T = atof(argv[++i]);
		} else if (strcmp(argv[i], "-ofn")==0) {
			params.outputFileName = argv[++i];
		} else if (strcmp(argv[i], "-cfn")==0) {
			params.cascadeFileName = argv[++i];
	} else if (strcmp(argv[i], "-mfn")==0) {
			params.modelFileName = argv[++i];
		} else if (strcmp(argv[i], "-wl")==0) {
			params.writeLog = (atoi(argv[++i])==1)?true:false;
		} else if (strcmp(argv[i], "-finSp")==0) {
			params.finSp = (atoi(argv[++i])==1)?true:false;
		} else if (strcmp(argv[i], "-mu")==0) {
			mu_mean = atof(argv[++i]);
		} else if (strcmp(argv[i], "-alpha")==0) {
			alpha_mean = atof(argv[++i]);
		} else if (strcmp(argv[i], "-eta")==0) {
			eta_mean = atof(argv[++i]);
		} else if (strcmp(argv[i], "-beta")==0) {
			beta_mean = atof(argv[++i]);
		} else if (strcmp(argv[i], "-w_phi")==0) {
			model.omega_phi = atof(argv[++i]);
		} else if (strcmp(argv[i], "-w_kap")==0) {
			model.omega_kap = atof(argv[++i]);
		} else if (strcmp(argv[i], "-rnd")==0) {
			params.random_params = (atoi(argv[++i])==1)?true:false;
		} else if (strcmp(argv[i], "-sp")==0) {
			params.sparsity = atof(argv[++i]);
		}
	}


	model.mu.resize(params.node_count);
	model.alpha.resize(params.node_count);
	model.eta.resize(params.node_count);
	model.beta.resize(params.node_count);

	if (params.random_params)  {
		for (int i=0; i < params.node_count; i++) {
			model.mu[i] = rng.uniform(0,mu_mean*2);
			model.alpha[i] = rng.uniform(0,alpha_mean*2);
			model.eta[i] = rng.uniform(0,eta_mean*2);
			model.beta[i] = rng.uniform(0,beta_mean*2);
		}
	} else {
		for (int i = 0; i < params.node_count; i++) {
			model.mu[i] = mu_mean;
			model.alpha[i] = alpha_mean;
			model.eta[i] = eta_mean;
			model.beta[i] = beta_mean;
		}
	}
	params.mu_mean = mu_mean;

	network.setSize(params.node_count);
	// history.setSize(params.node_count);
	Simulator* sim = new Simulator();

	sim->simulate();
	LOGMSG("simulation finished");
	cout << "simulation finished" << newline;
	sim->analyze_cascades();
	cout << "cascade analysis finished!" << newline;
}
