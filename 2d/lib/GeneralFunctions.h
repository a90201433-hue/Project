#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H

#include <functional>
//#include "Limiters.h"
#include "Types.h"

State PhysicalFlux(const State& W, int dir);
State ExactFlux(const State& WL,
                const State& WR,
                int dir);
State ComputeNumericalFlux(const State& WL,
                           const State& WR,
                           int dir);
void ReconstructGodunov(const Field& W,
                        Field& W_L,
                        Field& W_R,
                        int dir);
State Minmod(const State& a,
             const State& b);
State Superbee(const State& a,
			   const State& b);
void ReconstructKolgan(const Field& W,
                       Field& W_L,
                       Field& W_R,
                       int dir);

void Reconstruct(const Field& W,
                 Field& W_L,
                 Field& W_R,
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
