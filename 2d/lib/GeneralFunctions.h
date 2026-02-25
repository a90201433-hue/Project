#ifndef _GENERAL_FUNCTIONS_H_
#define _GENERAL_FUNCTIONS_H

#include <functional>
//#include "Limiters.h"
#include "Types.h"



//2D, по индексам, только точный решатель
void SolveBoundProblem(Field W, 
		       		   Field& W_b,
		       		   Field& F,
               		   int dir);


void FindBoundValues(Field W, Field& W_b,
					 Field& F,
                     std::vector<double> x,
					 std::vector<double> y,
					 double dt,
                     int dir);


void Streams(Field W,
			 Field& F,
    		 int dir);

// закомментил блок WENO и вязкости
void GetFluxes(Field W,
			   Field& F,
			   std::vector<double> x,
    		   std::vector<double> y,
			   double dt, 
    		   int dir);

// void Euler(Field& W_new, const Field& W, 
// 		   LimiterFunction phi,
// 		   RecLimiterFunction func, 
// 		   int init_idx_x, 
// 		   int end_idx_x, 
// 		   int init_idx_y,
// 		   int end_idx_y, 
// 		   const std::vector<double>& x,
// 		   const std::vector<double>& y,
// 		   double dt);
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
