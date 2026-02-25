#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H_

#include "Types.h"

State PhysicalFlux(const State& W, int dir);
State ExactFlux(const State& WL,
                const State& WR,
                int dir);
State ComputeNumericalFlux(const State& WL,
                           const State& WR,
                           int dir);

void ComputeFluxes(const Field& W,
                   Field& F,
                   int dir);

void Euler(const Field& W,
		   Field& W_new,
		   const std::vector<double>& x, 
		   const std::vector<double>& y,
		   double dt);

void UpdateArrays(Field& W, 
				  Field W_new,
				  std::vector<double> x, 
				  std::vector<double> y,
				  double dt) ;
#endif
