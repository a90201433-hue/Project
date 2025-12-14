#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H

#include <functional>

using LimiterFunction = std::function<double(double)>;

typedef double (*RecFunc)(double a, double b);
double sgn(double x);
double minmod(double a, double b);
double vanLeer(double a, double b);
double superbee(double a, double b);

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
	std::string high_order_method,
	std::string solver,
	std::string TVD_solver,
	LimiterFunction phi,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, double dt, bool Viscous_flag
);

void RK3(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::string method,
	
	std::string solver,
	
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt, bool Viscous_flag
	
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
	double dt, bool Viscous_flag
);

void MacCORMACK(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,	
	std::vector<double> x, 
	bool Viscous_flag,
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

void Visc(
	std::vector<double> W_L, 
	std::vector<double> W_R,
	double x_L, 
	double x_R,
	double& q);


void UpdateArrays(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,
	std::vector<std::vector<double>> W_b,
	std::vector<std::vector<double>> F,
	std::string method,
	std::string high_order_method,
	std::string solver,
	std::string TVD_solver,
	LimiterFunction phi,
	bool Viscous_flag,
	std::string time_method,
	std::vector<double> x, double dt);

#endif
