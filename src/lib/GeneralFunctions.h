#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H

typedef double (*RecFunc)(double a, double b);
double minmod(double a, double b);


void Streams(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F
);

void LRStreams(
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
    	std::vector<std::vector<double>>& F_left,
    	std::vector<std::vector<double>>& F_right
);


void TimePredictor(
	std::vector<std::vector<double>>& W_new,
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	int init_idx,
	int end_idx,
	std::vector<double> x,
	double dt
);
void Euler(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
    	std::string method, 
	std::string solver,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, double dt
);

void RK3(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::string method,
	std::string solver,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt
);

void FindSlopes(std::vector<std::vector<double>> W, std::vector<std::vector<double>>& Slope);

void ReconstructValues(
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>> Slope,
	std::vector<std::vector<double>>& W_L,
	std::vector<std::vector<double>>& W_R, 
	RecFunc function
);

void SolveBoundProblem(std::vector<std::vector<double>> W, 
		       std::vector<std::vector<double>>& W_b,
		       std::vector<std::vector<double>>& F,
		       std::string solver);

void SolveBoundProblem(std::vector<std::vector<double>> W_L, 
		       std::vector<std::vector<double>> W_R,
		       std::vector<std::vector<double>>& W_b,
		       std::vector<std::vector<double>>& F,
		       std::string solver);

void FindBoundValues(
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>>& W_b,
	std::vector<std::vector<double>>& F,
	std::vector<double> x, 
	double dt, 
	std::string method,
	std::string solver
);

void GetFluxes(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F,
	std::string method,
	std::string solver,
	std::vector<double> x,
	double dt
);

void MacCORMACK(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,	
	std::vector<double> x, 
	double dt, 
	std::string solver);

void DOperator(
	std::vector<std::vector<double>>& DU,
	std::vector<std::vector<double>> U,
	double Q);

void NOperator(
	std::vector<std::vector<double>>& NU, 
	std::vector<std::vector<double>> U,
	double Q);

void UpdateArrays(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,
	std::vector<std::vector<double>> W_b,
	std::vector<std::vector<double>> F,
	std::string method,
	std::string solver,
	std::string time_method,
	std::vector<double> x,double dt, double mu0);

#endif
