#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H

#include <functional>
//#include "Limiters.h"
#include "Types.h"




//2D, по индексам, только точный решатель
void SolveBoundProblem(Field W, 
		       Field& W_b,
		       Field& F,
		       std::string solver,
               int dir);
void FindBoundValues(Field W, 
					 Field& W_b,
					 Field& F,
                     std::vector<double> x,
					 std::vector<double> y,
					 double dt,
					 std::string method,
					 std::string solver,
					 RecLimiterFunction func,
                     int dir);
void Streams(
	Field W,
	Field& F,
    int dir);

void GetFluxes(// закомментил блок WENO и вязкости
	Field W,
	Field& F,
	std::string method,
	std::string solver,
	RecLimiterFunction func,
	std::vector<double> x,
    std::vector<double> y,
	double dt, 
	bool Viscous_flag,
    int dir);

void Euler(
	Field& W_new, 
	Field W, 
    	std::string method, 
	std::string high_order_method,
	std::string solver,
	std::string TVD_solver,
	LimiterFunction phi,
	RecLimiterFunction func,
	int init_idx_x, 
	int end_idx_x, 
    int init_idx_y, 
	int end_idx_y, 
	std::vector<double> x, 
    std::vector<double> y, 
	double dt, bool Viscous_flag
    );

void UpdateArrays(Field& W, Field W_new,
				  std::string method,
				  std::string high_order_method,
				  std::string solver,
				  std::string TVD_solver,
				  LimiterFunction phi,
				  RecLimiterFunction func,
				  bool Viscous_flag,
				  std::string time_method,
				  std::vector<double> x, 
				  std::vector<double> y,
				  double dt);
#endif
