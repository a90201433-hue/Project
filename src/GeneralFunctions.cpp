#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "lib/RiemannSolver.h"
#include "lib/BoundCond.h"
#include "lib/TransformValues.h"
#include "lib/Init.h"
#include "lib/GeneralFunctions.h"

extern double gamm, L;
extern int N, fict;

typedef double (*RecFunc)(double a, double b);

double minmod(double a, double b) {
	if (a * b <= 0) {
		return 0.0;
	}
	else if (a * a <= a * b) {
		return a;
	}
	else if (b * b < a * b) {
		return b;
	}
	return 0.0;
}

typedef void (*FluxFunc) (
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	std::vector<std::vector<double>>& F_left,
	std::vector<std::vector<double>>& F_right);

void GodunovStreams(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F) {

	InitialZeros(F, 3);
	for (int i = fict; i < N + fict; i++) {
		F[i][0] = W[i][0] * W[i][1]; // rho*u

		F[i][1] = W[i][0] * std::pow(W[i][1], 2) + W[i][2]; // rho*u^2 + P
		
		double E = 0.5 * std::pow(W[i][1], 2) + W[i][2] / (W[i][0] * (gamm - 1.0));
		F[i][2] = W[i][1] * (W[i][0] * E + W[i][2]); // u(rho*E + P)
	}
}

void ConservativeFlux(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
    	std::vector<std::vector<double>>& F_left,
    	std::vector<std::vector<double>>& F_right) {

	GodunovStreams(W, F_left);
    	F_right = F_left; 
}

void NonConservativeFlux(
    	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
    	std::vector<std::vector<double>>& F_left,
    	std::vector<std::vector<double>>& F_right) {

    	GodunovStreams(W_L, F_right);  
    	GodunovStreams(W_R, F_left);   
}

typedef void (*TimeFunc)(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	FluxFunc ComputeFlux,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt);

void Euler(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
    	std::string method, 
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt) {

	std::vector<std::vector<double>> F(N + 2 * fict);
	InitialZeros(F, 3);
    
	GetFluxes(W, F, method, x, dt);
	for (int i = init_idx; i < end_idx; i++) {
		double dx = x[i + 1] - x[i];

		W_new[i][0] = std::max(1e-7, W[i][0] 
			    - dt/dx * (F[i + 1][0] - F[i][0]));

		W_new[i][1] = (1.0 / W_new[i][0]) * (W[i][0] * W[i][1] 
			    -  dt/dx * (F[i + 1][1] - F[i][1]));	

		W_new[i][2] = std::max(1e-7, (gamm - 1.0) 
			    * (0.5 * (W[i][0] * std::pow(W[i][1], 2) 
			    - W_new[i][0] * std::pow(W_new[i][1], 2)) 
			    + W[i][2] / (gamm - 1.0) - dt / dx * (F[i + 1][2] - F[i][2])));	
	}
}

void RK3(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::string method,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt) {

	std::vector<std::vector<double>> W1(N + 2 * fict - 1);
	InitialZeros(W1, 3);
	std::vector<std::vector<double>> W2(N + 2 * fict - 1);
	InitialZeros(W2, 3);
	std::vector<std::vector<double>> F(N + 2 * fict);
	InitialZeros(F, 3);
	

	GetFluxes(W, F, method, x, dt);
	for (int i = init_idx; i < end_idx; i++) {
		double dx = x[i + 1] - x[i];

		W1[i][0] = std::max(1e-7, W[i][0] 
			 - dt/dx * (F[i + 1][0] - F[i][0]));

		W1[i][1] = (1.0 / W1[i][0]) * (W[i][0] * W[i][1] 
			 - dt / dx * (F[i + 1][1] - F[i][1]));

		W1[i][2] = std::max(1e-7, (gamm - 1.0) * (0.5 * (W[i][0] * std::pow(W[i][1], 2) 
			 - W1[i][0] * std::pow(W1[i][1], 2))
			 + W[i][2] / (gamm - 1.0) 
			 - dt / dx * (F[i + 1][2] - F[i][2])));
	}
	BoundCond(W1);

	GetFluxes(W1, F, method, x, dt);
	for (int i = init_idx; i < end_idx; i++) {
		double dx = x[i + 1] - x[i];

		double W1_temp_0 = std::max(1e-7, W1[i][0] 
				  - dt/dx * (F[i + 1][0] - F[i][0]));
        	W2[i][0] = std::max(1e-7, 0.75 * W[i][0] + 0.25 * W1_temp_0);


		double W1_temp_1 = (1.0 / W1_temp_0) * (W1[i][0] * W1[i][1] 
				  - dt/dx * (F[i + 1][1] - F[i][1]));
        	W2[i][1] = 0.75 * W[i][1] + 0.25 * W1_temp_1;

		double W1_temp_2 = std::max(1e-7, (gamm - 1.0) * (0.5 * (W1[i][0] * std::pow(W1[i][1], 2) 
				  - W1_temp_0 * std::pow(W1_temp_1, 2))
				  + W1[i][2] / (gamm - 1.0) 
				  - dt/dx * (F[i + 1][2] - F[i][2])));
		W2[i][2] = std::max(1e-7, 0.75 * W[i][2] + 0.25 * W1_temp_2);
	}
	BoundCond(W2);

	GetFluxes(W2, F, method, x, dt);
	for (int i = init_idx; i < end_idx; i++) {
		double dx = x[i + 1] - x[i];
			
		double W2_temp_0 = std::max(1e-7, W2[i][0] 
				  - dt/dx * (F[i + 1][0] - F[i][0]));
        	W_new[i][0] = std::max(1e-7, (1.0/3.0) * W[i][0] + (2.0/3.0) * W2_temp_0);

		double W2_temp_1 = (1.0 / W2_temp_0) * (W2[i][0] * W2[i][1] 
				  - dt/dx * (F[i + 1][1] - F[i][1]));
        	W_new[i][1] = (1.0/3.0) * W[i][1] + (2.0/3.0) * W2_temp_1;

		double W2_temp_2 = std::max(1e-7, (gamm - 1.0) * (0.5 * (W2[i][0] * std::pow(W2[i][1], 2) 
				  - W2_temp_0 * std::pow(W2_temp_1, 2))
				  + W2[i][2] / (gamm - 1.0)
				   - dt/dx * (F[i + 1][2] - F[i][2])));
		W_new[i][2] = std::max(1e-7, (1.0/3.0) * W[i][2] + (2.0/3.0) * W2_temp_2);
	}
}


void TimeRodionovPredictor(
	std::vector<std::vector<double>>& W_new, 
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>> W_L,
	std::vector<std::vector<double>> W_R,
	FluxFunc ComputeFlux, 
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt) {

	std::vector<std::vector<double>> F_left(N + 2 * fict); 
	InitialZeros(F_left, 3);
	std::vector<std::vector<double>> F_right(N + 2 * fict); 
	InitialZeros(F_right, 3);
    
	ComputeFlux(W, W_L, W_R, F_left, F_right); 
	//W_tilde 
	for (int i = init_idx; i < end_idx; i++) {
		double dx = x[i + 1] - x[i];
		W_new[i][0] = std::max(1e-7, W[i][0] 
			    - dt/dx * (F_right[i + 1][0] - F_left[i][0]));

		W_new[i][1] = (1.0 / W_new[i][0]) * (W[i][0] * W[i][1] 
			    - dt / dx * (F_right[i + 1][1] - F_left[i][1]));

		W_new[i][2] = std::max(1e-7, (gamm - 1.0) * (0.5 * (W[i][0] 
			    * std::pow(W[i][1], 2) 
			    - W_new[i][0] * std::pow(W_new[i][1], 2)) 
			    + W[i][2] / (gamm - 1.0) 
			    - dt / dx * (F_right[i + 1][2] - F_left[i][2])));
	}
}


void FindSlopes(std::vector<std::vector<double>> W, std::vector<std::vector<double>>& Slope){
	InitialZeros(Slope, 3);
	for (int j = 0; j < 3; j++) {
		for (int i = 1; i < N + 2 * fict - 1; i++) {
			Slope[i][j] = W[i][j] - W[i - 1][j];
		}
	}
}

void ReconstructValues(
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>> Slope,
	std::vector<std::vector<double>>& W_L,
	std::vector<std::vector<double>>& W_R, 
	RecFunc function
) {
	for (int j = 0; j < 3; j++){
		for (int i = fict; i < N + fict; i++) {
			W_R[i][j] = W[i][j] - 0.5 * function(Slope[i][j], Slope[i + 1][j]);
		
			W_L[i][j] = W[i - 1][j] + 0.5 * function(Slope[i - 1][j], Slope[i][j]);
		}
	}
}

double ENOPolinom(double f1, double f2, double f3, double x) {
	double h = L / (N - 1);
	double result = ((11 * f1 - 7 * f2 + 2 * f3) * std::pow(h, 2) 
		- 12 * (f1 - 1.5 * f2 + 0.5 * f3) * x * h
		+ 3 * std::pow(x, 2) * (f1 - 2 * f2 + f3)) / (6.0 * std::pow(h, 2));
	return result;
}	

void ENO(
	std::vector<std::vector<double>> U,
	std::vector<std::vector<double>> Slope,
	std::vector<std::vector<double>>& U_L,
	std::vector<std::vector<double>>& U_R) {	

	double h = L / (N - 1);

	for (int j = 0; j < 3; j++) {
		for (int i = fict - 1; i < N + fict - 1; i++) {
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
	
		for (int i = fict; i < N + fict; i++) {
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
	std::vector<std::vector<double>> U,
	std::vector<std::vector<double>>& U_L,
	std::vector<std::vector<double>>& U_R) {	

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
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F) {	

	std::vector<std::vector<double>> W_L(N + 2 * fict);
	std::vector<std::vector<double>> W_R(N + 2 * fict);
	InitialZeros(W_L, 3);
	InitialZeros(W_R, 3);
	
	std::vector<std::vector<double>> U_L(N + 2 * fict);
	std::vector<std::vector<double>> U_R(N + 2 * fict);
	InitialZeros(U_L, 3);
	InitialZeros(U_R, 3);

	//Reconstuct WENO
	WENO(W, W_L, W_R);	
		
	ConvertWtoU(W_L, U_L);
	ConvertWtoU(W_R, U_R);
	
	std::vector<std::vector<double>> F_L(N + 2 * fict);
	std::vector<std::vector<double>> F_R(N + 2 * fict);	

	InitialZeros(F_L, 3);
	InitialZeros(F_R, 3);

	//Godunov streams
	GodunovStreams(W_L, F_L); //F(U_L)
	GodunovStreams(W_R, F_R); //F(U_R)
	
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

void FindBoundValues(
	std::vector<std::vector<double>> W, 
	std::vector<std::vector<double>>& W_b, 
	std::vector<double> x, 
	double dt, 
	std::string method
){
	std::vector<std::vector<double>> Slope(N + 2 * fict);
	InitialZeros(Slope, 3);
	std::vector<std::vector<double>> W_L(N + 2 * fict);
	InitialZeros(W_L, 3);
	std::vector<std::vector<double>> W_R(N + 2 * fict);
	InitialZeros(W_R, 3);
	std::vector<std::vector<double>> U(N + 2 * fict - 1);
	InitialZeros(U, 3);
	std::vector<std::vector<double>> U_L(N + 2 * fict);
	InitialZeros(U_L, 3);
	std::vector<std::vector<double>> U_R(N + 2 * fict);
	InitialZeros(U_R, 3);
	std::vector<std::vector<double>> W_tilde(N + 2 * fict - 1);
	InitialZeros(W_tilde, 3);
	std::vector<std::vector<double>> W_half(N + 2 * fict - 1);
	InitialZeros(W_half, 3);
	std::vector<std::vector<double>> U_half(N + 2 * fict - 1);
	InitialZeros(U_half, 3);

	if (method == "Godunov"){

		// 1. Ищем звездные значения
		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W[i - 1], W[i], W_b[i], 1e-6);
		}
		// 2. Ищем значение на границе из задаче о распаде разрыва
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W[i - 1], W[i], W_b[i], 0.0, 1.0);
		}
		return;
	}

	else if (method == "Kolgan"){
	
		FindSlopes(W, Slope);	
		ReconstructValues(W, Slope, W_L, W_R, &minmod);

		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6);
		}
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);
		}
		return;		
	}
	
	else if (method == "Kolgan2"){

		ConvertWtoU(W, U);
		FindSlopes(U, Slope);
		ReconstructValues(U, Slope, U_L, U_R, &minmod);
		ConvertUtoW(W_L, U_L);
		ConvertUtoW(W_R, U_R);

		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6);	
		}
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);	
		}
		return;		
	}

	else if (method == "Rodionov"){
		
	FindSlopes(W, Slope);
		ReconstructValues(W, Slope, W_L, W_R, &minmod);
		TimeRodionovPredictor(W_tilde, W, W_L, W_R, &NonConservativeFlux, fict, N + fict - 1, x, dt);
		BoundCond(W_tilde);
		
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < N + 2 * fict - 1; i++){
				W_half[i][j] = 0.5 * (W[i][j] + W_tilde[i][j]);	
			}
		}
		
		ReconstructValues(W_half, Slope, W_L, W_R, &minmod);
		
		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6);
		}
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);
		}
		return;	
	}
	
	else if (method == "Rodionov2"){
		
		ConvertWtoU(W, U);
		FindSlopes(U, Slope);	
		ReconstructValues(U, Slope, U_L, U_R, &minmod);
		ConvertUtoW(W_L, U_L);
		ConvertUtoW(W_R, U_R);
		TimeRodionovPredictor(W_tilde, W, W_L, W_R, &NonConservativeFlux, fict, N + fict - 1, x, dt);
		BoundCond(W_tilde);
		
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < N + 2 * fict - 1; i++){
				W_half[i][j] = 0.5 * (W[i][j] + W_tilde[i][j]);
			}
		}
		
		ConvertWtoU(W_half, U_half);
		ReconstructValues(U_half, Slope, U_L, U_R, &minmod);
		ConvertUtoW(W_L, U_L);
		ConvertUtoW(W_R, U_R);

		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6);
		}
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);
		}
		return;	
	}

	else if (method == "ENO") {

		
		FindSlopes(W, Slope);
		ENO(W, Slope, W_L, W_R);

		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6);	
		}
		for (int i = fict; i < N + fict; i++) {
			
			W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);	
		}
		return;
	}
}
	
void GetFluxes(
	std::vector<std::vector<double>> W,
	std::vector<std::vector<double>>& F,
	std::string method,
	std::vector<double> x,
	double dt
) {
	if (method == "WENO") {
		WENOStreams(W, F);
	} else {
		std::vector<std::vector<double>> W_b(N + 2 * fict);
		InitialZeros(W_b, 3);
		FindBoundValues(W, W_b, x, dt, method); 
		GodunovStreams(W_b, F);
	}
}
void UpdateArrays(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,
	std::vector<std::vector<double>> W_b,
	std::vector<std::vector<double>> F,
	std::vector<double> x, 
	double dt, 
	std::string method,
	std::string time_method
){
	std::vector<std::vector<double>> W_L(N + 2 * fict);
	std::vector<std::vector<double>> W_R(N + 2 * fict);

	if (time_method == "RK3") {
		RK3(W_new, W, method, fict, N + fict - 1, x, dt);
	} else {
		Euler(W_new, W, method, fict, N + fict - 1, x, dt);
	}
	//std::cout << "Time method = " << time_method << std::endl;

	for (int i = fict; i < N + fict - 1; i++) {
		W[i] = W_new[i];
	}
}
