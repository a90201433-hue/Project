#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "lib/ENO.h"
#include "lib/TransformValues.h"
#include "lib/GeneralFunctions.h"
#include "lib/Init.h"

extern double gamm, L;
extern int N, fict;

using DataArray = std::vector<std::vector<double>>;

std::vector<double> WENOWeight(double f1, double f2, double f3, double f4, double f5) {
	double beta0 = 13 * std::pow(f3 - 2 * f4 + f5, 2)/12 + std::pow(3 * f3 - 4 * f4 + f5, 2)/4; //i, i+1, i+2
	double beta1 = 13 * std::pow(f2 - 2 * f3 + f4, 2)/12 + std::pow(f2 - f4, 2)/4; //i-1, i, i+1
	double beta2 = 13 * std::pow(f1 - 2 * f2 + f3, 2)/12 + std::pow(f1 - 4 * f2 + 3 * f3, 2)/4; //i-2, i-1, i
	 
	
	double d0 = 0.3;
	double d1 = 0.6;
	double d2 = 0.1;
	
	double eps = 1e-7;

	double alpha0 = d0 / std::pow(eps + beta0, 2);
	double alpha1 = d1 / std::pow(eps + beta1, 2);
	double alpha2 = d2 / std::pow(eps + beta2, 2);
	
	std::vector<double> all_weights = {0.0, 0.0, 0.0};
	all_weights = {alpha0/(alpha0 + alpha1 + alpha2), 
		       alpha1/(alpha0 + alpha1 + alpha2), 
		       alpha2/(alpha0 + alpha1 + alpha2)};
	for (int i = 0; i < 3; i++){
		if (all_weights[i] < 1e-6)
			all_weights[i] = 0.0;
	}	
	return all_weights;
}

void WENO(
	DataArray U,
	DataArray& U_L,
	DataArray& U_R) {	

	double value1, value2, value3;
	double h = L / (N - 1);
	std::vector<double> coeffs = {0.0, 0.0, 0.0};
	
	for (int j = 0; j < 3; j++) {
		for (int i = fict - 1; i < N + fict - 1; i++) {
			
			coeffs = WENOWeight(U[i - 2][j], U[i - 1][j], U[i][j], U[i + 1][j], U[i + 2][j]);
			

			if (j == 0 || j == 2) { 
				value1 = ENOPolinom(U[i][j], U[i + 1][j], U[i + 2][j], h);
				value2 = ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], 2 * h);
				value3 = ENOPolinom(U[i - 2][j], U[i - 1][j], U[i][j], 3 * h);
				

				U_L[i + 1][j] = coeffs[0] * value1 + coeffs[1] * value2 + coeffs[2] * value3;
				U_L[i + 1][j] = std::max(1e-6, U_L[i + 1][j]);

			} else {
					
				U_L[i + 1][j] = coeffs[0] * ENOPolinom(U[i][j], U[i + 1][j], U[i + 2][j], h) 
				      	+ coeffs[1] * ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], 2 * h)
				      	+ coeffs[2] * ENOPolinom(U[i - 2][j], U[i - 1][j], U[i][j], 3 * h);
				
			}
			

		}
		
		for (int i = fict; i < N + fict; i++) {
			coeffs = WENOWeight(U[i - 2][j], U[i - 1][j], U[i][j], U[i + 1][j], U[i + 2][j]);
			if (j == 0 || j == 2) { 
				value1 = ENOPolinom(U[i][j], U[i + 1][j], U[i + 2][j], 0);
				value2 = ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], h);
				value3 = ENOPolinom(U[i - 2][j], U[i - 1][j], U[i][j], 2 * h);
				U_R[i][j] = coeffs[0] * value1 + coeffs[1] * value2 + coeffs[2] * value3;
				U_R[i][j] = std::max(1e-6, U_R[i][j]);

			} else {
				U_R[i][j] = coeffs[0] * ENOPolinom(U[i][j], U[i + 1][j], U[i + 2][j], 0) 
				  	+ coeffs[1] * ENOPolinom(U[i - 1][j], U[i][j], U[i + 1][j], h)
				  	+ coeffs[2] * ENOPolinom(U[i - 2][j], U[i - 1][j], U[i][j], 2 * h);
			}
			
		}
	}
}

void WENOStreams(
	DataArray W,
	DataArray& F) {	

	DataArray W_L(N + 2 * fict);
	DataArray W_R(N + 2 * fict);
	InitialZeros(W_L, 3);
	InitialZeros(W_R, 3);
	
	DataArray U_L(N + 2 * fict);
	DataArray U_R(N + 2 * fict);
	InitialZeros(U_L, 3);
	InitialZeros(U_R, 3);

	//Reconstuct WENO
	WENO(W, W_L, W_R);	
		
	ConvertWtoU(W_L, U_L);
	ConvertWtoU(W_R, U_R);
	
	DataArray F_L(N + 2 * fict);
	DataArray F_R(N + 2 * fict);	

	InitialZeros(F_L, 3);
	InitialZeros(F_R, 3);

	//Godunov streams
	Streams(W_L, F_L); //F(U_L)
	Streams(W_R, F_R); //F(U_R)
	
	std::vector<double> c_L(N + 2 * fict);
	std::vector<double> c_R(N + 2 * fict);
	for (int i = 0; i < N + 2 * fict; i++) {
		c_L[i] = 0.0;
		c_R[i] = 0.0;
	}
	for (int i = fict; i < N + fict; i++) {
		c_L[i] = std::sqrt(gamm * W_L[i][2]/W_L[i][0]);
		c_R[i] = std::sqrt(gamm * W_R[i][2]/W_R[i][0]);
	}

	std::vector<double> alpha(N + 2 * fict);
	for (int i = 0; i < N + 2 * fict; i++) {
		alpha[i] = 0.0;
	}
	for (int i = 0; i < N + 2 * fict; i++) {
		alpha[i] = std::max(std::abs(W_L[i][1]) + c_L[i], std::abs(W_R[i][1]) + c_R[i]);
		//alpha[i] = std::max(std::abs(c_R[i]), std::abs(c_L[i]));
	}

	//WENO streams
	for (int i = fict; i < N + fict; i++) {
		for (int j = 0; j < 3; j++) {
			F[i][j] = 0.5 * ((F_R[i][j] + F_L[i][j]) - alpha[i] * (U_R[i][j] - U_L[i][j]));
		}
	}
}
