#include <vector>
#include <cmath>
#include <iostream>

extern double gamm, Lx;
extern int Nx, fict;

using DataArray = std::vector<std::vector<double>>;


double ENOPolinom(double f1, double f2, double f3, double x) {
	double h = Lx / (Nx - 1);
	double result = ((11 * f1 - 7 * f2 + 2 * f3) * std::pow(h, 2) 
		- 12 * (f1 - 1.5 * f2 + 0.5 * f3) * x * h
		+ 3 * std::pow(x, 2) * (f1 - 2 * f2 + f3)) / (6.0 * std::pow(h, 2));
	return result;
}

void ENO(
	DataArray U,
	DataArray Slope,
	DataArray& U_L,
	DataArray& U_R) {	

	double h = Lx / (Nx - 1);

	for (int j = 0; j < 3; j++) {
		for (int i = fict - 1; i < Nx + fict - 1; i++) {
			if (std::abs(Slope[i][j]) < std::abs(Slope[i + 1][j])) { 
				if (std::abs(Slope[i][j] - Slope[i - 1][j]) < std::abs(Slope[i + 1][j] - Slope[i][j])) {	
					U_L[i + 1][j] = ENOPolinom(U[i - 2][j], U[i - 1][j], U[i][j], 3 * h); 
				} else {
					U_L[i + 1][j] = ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], 2 * h);
				}										
			} else {
				if (std::abs(Slope[i + 2][j] - Slope[i + 1][j]) < std::abs(Slope[i + 1][j] - Slope[i][j])) { 	
					U_L[i + 1][j] = ENOPolinom(U[i][j], U[i + 1][j], U[i + 2][j], h);
							
				} else {
					U_L[i + 1][j] = ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], 2 * h);				
				}								
			}
			if (j == 0 || j == 2) {
				U_L[i + 1][j] = std::max(1e-6, U_L[i + 1][j]);
			}
						
		}
	
		for (int i = fict; i < Nx + fict; i++) {
			if (std::abs(Slope[i][j]) < std::abs(Slope[i + 1][j])) { 
				if (std::abs(Slope[i][j] - Slope[i - 1][j]) < std::abs(Slope[i + 1][j] - Slope[i][j])) {	
					U_R[i][j] = ENOPolinom(U[i - 2][j], U[i - 1][j], U[i][j], 2 * h);
				} else {			
					U_R[i][j] = ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], h);		
				}										
			} else {
				if (std::abs(Slope[i + 2][j] - Slope[i + 1][j]) < std::abs(Slope[i + 1][j] - Slope[i][j])) { 
					U_R[i][j] = ENOPolinom(U[i][j], U[i + 1][j], U[i + 2][j], 0);					
				} else {
					U_R[i][j] = ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], h);	
				}							
			}	
			if (j == 0 || j == 2) {
				U_R[i][j] = std::max(1e-6, U_R[i][j]);
			}		
		}
	}
}
