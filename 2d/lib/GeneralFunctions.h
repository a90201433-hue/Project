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


void FindBoundValues(Field W, 
					 Field& W_b,
					 Field& F,
                     std::vector<double> x,
					 std::vector<double> y,
					 double dt,
					 RecLimiterFunction func,
                     int dir);


void Streams(Field W,
			 Field& F,
    		 int dir);

// закомментил блок WENO и вязкости
void GetFluxes(Field W,
			   Field& F,
			   RecLimiterFunction func,
			   std::vector<double> x,
    		   std::vector<double> y,
			   double dt, 
    		   int dir);

void Euler(Field& W_new, 
		   Field W, 
		   LimiterFunction phi,
		   RecLimiterFunction func,
		   int init_idx_x, 
		   int end_idx_x, 
    	   int init_idx_y, 
		   int end_idx_y, 
		   std::vector<double> x, 
    	   std::vector<double> y, 
		   double dt);

void UpdateArrays(Field& W, 
				  Field W_new,
				  LimiterFunction phi,
				  RecLimiterFunction func,
				  std::vector<double> x, 
				  std::vector<double> y,
				  double dt);
#endif
