#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H

typedef double (*RecFunc)(double a, double b);
double minmod(double a, double b);


typedef void (*FluxFunc) (
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	std::vector<std::vector<double>>& F_left,
	std::vector<std::vector<double>>& F_right
);
void GodunovStreams(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F
);
void ConservativeFlux(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
    	std::vector<std::vector<double>>& F_left,
    	std::vector<std::vector<double>>& F_right
);
void NonConservativeFlux(
    	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
    	std::vector<std::vector<double>>& F_left,
    	std::vector<std::vector<double>>& F_right
);


typedef void (*TimeFunc)(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	FluxFunc ComputeFlux,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt
);
void TimeRodionovPredictor(
	std::vector<std::vector<double>>& W_new,
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	FluxFunc ComputeFlux,
	int init_idx,
	int end_idx,
	std::vector<double> x,
	double dt
);
void Euler(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
    	std::string method, 
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt
);
void RK3(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::string method,
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


double ENOPolinom(double f1, double f2, double f3, double x);
void ENO(
	std::vector<std::vector<double>> U,
	std::vector<std::vector<double>> Slope,
	std::vector<std::vector<double>>& U_L,
	std::vector<std::vector<double>>& U_R
);
std::vector<double> WENOWeight(double f1, double f2, double f3, double f4, double f5);
void WENO(
	std::vector<std::vector<double>> U,
	std::vector<std::vector<double>>& U_L,
	std::vector<std::vector<double>>& U_R
);
void WENOStreams(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F
);


void FindBoundValues(
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>>& W_b, 
	std::vector<double> x, 
	double dt, 
	std::string method
);

void GetFluxes(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F,
	std::string method,
	std::vector<double> x,
	double dt
);

void MacCORMACK(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,	
	std::vector<double> x, 
	double dt);

void UpdateArrays(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,
	std::vector<std::vector<double>> W_b,
	std::vector<std::vector<double>> F,
	std::vector<double> x, 
	double dt, 
	std::string method,
	std::string time_method
);

#endif
