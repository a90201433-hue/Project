#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "lib/RiemannSolver.h"
#include "lib/ApproxSolvers.h"
#include "lib/BoundCond.h"
#include "lib/TransformValues.h"
#include "lib/Init.h"
#include "lib/GeneralFunctions.h"
#include "lib/ENO.h"
#include "lib/WENO.h"
#include "lib/Limiters.h"

extern double gamm, L;
extern int N, fict;

using DataArray = std::vector<std::vector<double>>;

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

void Streams(
	DataArray W,
	DataArray& F) {

	InitialZeros(F, 3);
	for (int i = fict; i < N + fict; i++) {
		F[i][0] = W[i][0] * W[i][1]; // rho*u

		F[i][1] = W[i][0] * std::pow(W[i][1], 2) + W[i][2]; // rho*u^2 + P
		
		double E = 0.5 * std::pow(W[i][1], 2) + W[i][2] / (W[i][0] * (gamm - 1.0));
		F[i][2] = W[i][1] * (W[i][0] * E + W[i][2]); // u(rho*E + P)
	}
}


void LRStreams(
	DataArray W_L,
	DataArray W_R,
    	DataArray& F_left,
    	DataArray& F_right) {

    	Streams(W_L, F_right);  
    	Streams(W_R, F_left);   
}


void Euler(
	DataArray& W_new, 
	DataArray W, 
    std::string method, 
	std::string high_order_method,
	std::string solver,
	std::string TVD_solver,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt) {

	DataArray F(N + 2 * fict);
	InitialZeros(F, 3);
    
	if (method == "TVD") {
		DataArray F_low(N + 2 * fict, {0.0, 0.0, 0.0});
		GetFluxes(W, F_low, "Godunov", TVD_solver, x, dt);

		DataArray F_high(N + 2 * fict, {0.0, 0.0, 0.0});
		GetFluxes(W, F_high, high_order_method, TVD_solver, x, dt);

		std::vector<double> r_array(3, 0.0);

		for (int i = fict; i < N + fict; i++) {
			r_array = FindR(W[i - 1], W[i], W[i+1]);
			for (int j = 0; j < 3; j++) {
				F[i][j] = F_low[i][j] - CHARM(r_array[j])*(F_low[i][j] - F_high[i][j]);
			}
		}
	}
	else {
		GetFluxes(W, F, method, solver, x, dt);
	}

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
	DataArray& W_new, 
	DataArray W, 
	std::string method,
	std::string solver,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt) {

	DataArray W1(N + 2 * fict - 1);
	InitialZeros(W1, 3);
	DataArray W2(N + 2 * fict - 1);
	InitialZeros(W2, 3);
	DataArray F(N + 2 * fict);
	InitialZeros(F, 3);
	

	GetFluxes(W, F, method, solver, x, dt);
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

	GetFluxes(W1, F, method, solver, x, dt);
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

	GetFluxes(W2, F, method,solver, x, dt);
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


void TimePredictor(
	DataArray& W_new, 
	DataArray W, 
	DataArray W_L,
	DataArray W_R,
	int init_idx, 
	int end_idx, 
	std::vector<double> x, 
	double dt) {

	DataArray F_left(N + 2 * fict); 
	InitialZeros(F_left, 3);
	DataArray F_right(N + 2 * fict); 
	InitialZeros(F_right, 3);
    
	LRStreams(W_L, W_R, F_left, F_right); 
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


void FindSlopes(std::vector<std::vector<double>> W, 
				std::vector<std::vector<double>>& Slope) {
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


void SolveBoundProblem(std::vector<std::vector<double>> W, 
		       std::vector<std::vector<double>>& W_b,
		       std::vector<std::vector<double>>& F,
		       std::string solver){
	
	if (solver == "Exact") {

		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W[i - 1], W[i], W_b[i], 1e-6);
		}
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W[i - 1], W[i], W_b[i], 0.0, 1.0);
		}	
		return;
	}

	if (solver == "HLL") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = HLL(W[i - 1], W[i]);	
		}
		return;
	}

	if (solver == "HLLC") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = HLLC(W[i - 1], W[i]);	
		}
		return;
	}

	if (solver == "Rusanov") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = Rusanov(W[i - 1], W[i]);	
		}
		return;
	}

	if (solver == "Roe") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = Roe(W[i - 1], W[i]);	
		}
		return;
	}

	if (solver == "Osher") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = Osher(W[i - 1], W[i]);	
		}
		return;
	}

}

void SolveBoundProblem(std::vector<std::vector<double>> W_L, 
		       		   std::vector<std::vector<double>> W_R,
		       		   std::vector<std::vector<double>>& W_b,
		       		   std::vector<std::vector<double>>& F,
		       		   std::string solver){
	
	if (solver == "Exact") {

		for (int i = fict; i < N + fict; i++) {
			NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6);
		}
		for (int i = fict; i < N + fict; i++) {
			W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);
		}	
		return;
	}

	if (solver == "HLL") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = HLL(W_L[i], W_R[i]);	
		}
		return;
	}

	if (solver == "HLLC") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = HLLC(W_L[i], W_R[i]);	
		}
		return;
	}

	if (solver == "Rusanov") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = Rusanov(W_L[i], W_R[i]);	
		}
		return;
	}

	if (solver == "Roe") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = Roe(W_L[i], W_R[i]);	
		}
		return;
	}

	if (solver == "Osher") {
		for (int i = fict; i < N + fict; i++) {
			F[i] = Osher(W_L[i], W_R[i]);	
		}
		return;
	}


}


void FindBoundValues(std::vector<std::vector<double>> W, 
					 std::vector<std::vector<double>>& W_b,
					 std::vector<std::vector<double>>& F,
					 std::vector<double> x, 
					 double dt, 
					 std::string method,
					 std::string solver) {

	std::vector<std::vector<double>> Slope(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_L(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_R(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_L(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_R(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_tilde(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_half(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_half(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	

	if (method == "Godunov") {
		SolveBoundProblem(W, W_b, F, solver);
		return;
	}

	else if (method == "Kolgan") {
	
		FindSlopes(W, Slope);	
		ReconstructValues(W, Slope, W_L, W_R, &minmod);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);
		return;		
	}
	
	else if (method == "Kolgan2") {

		ConvertWtoU(W, U);
		FindSlopes(U, Slope);
		ReconstructValues(U, Slope, U_L, U_R, &minmod);
		ConvertUtoW(W_L, U_L);
		ConvertUtoW(W_R, U_R);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);

		return;		
	}

	else if (method == "Rodionov") {
		
		FindSlopes(W, Slope);
		ReconstructValues(W, Slope, W_L, W_R, &minmod);
		TimePredictor(W_tilde, W, W_L, W_R, fict, N + fict - 1, x, dt);
		BoundCond(W_tilde);
		
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < N + 2 * fict - 1; i++){
				W_half[i][j] = 0.5 * (W[i][j] + W_tilde[i][j]);	
			}
		}
		
		ReconstructValues(W_half, Slope, W_L, W_R, &minmod);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);
		return;	
	}
	
	else if (method == "Rodionov2") {
		
		ConvertWtoU(W, U);
		FindSlopes(U, Slope);	
		ReconstructValues(U, Slope, U_L, U_R, &minmod);
		ConvertUtoW(W_L, U_L);
		ConvertUtoW(W_R, U_R);
		TimePredictor(W_tilde, W, W_L, W_R, fict, N + fict - 1, x, dt);
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
		SolveBoundProblem(W_L, W_R, W_b, F, solver);
		return;	
	}

	else if (method == "ENO") {
		
		FindSlopes(W, Slope);
		ENO(W, Slope, W_L, W_R);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);
		return;
	}
}
	
void GetFluxes(
	DataArray W,
	DataArray& F,
	std::string method,
	std::string solver,
	std::vector<double> x,
	double dt
) {
	if (method == "WENO") {
		WENOStreams(W, F);
	} else {
		std::vector<std::vector<double>> W_b(N + 2 * fict);
		InitialZeros(W_b, 3);
		FindBoundValues(W, W_b, F, x, dt, method, solver);
		if (solver == "Exact"){
			Streams(W_b, F);
		}
	}
}

void MacCORMACK(
	DataArray& W,
	DataArray W_new,	
	std::vector<double> x, 
	bool Viscous_flag, double mu0,
	double dt,
	std::string solver){
	
	const size_t Central_num = N + 2 * fict - 1; 
	
	DataArray F(Central_num);
	InitialZeros(F, 3);
	Streams(W, F); 

	DataArray U(Central_num);
	ConvertWtoU(W, U);

	DataArray U_tilde(Central_num);
	InitialZeros(U_tilde, 3);
	

	for (int j = 0; j < 3; j++) {
		for (int i = fict; i < N + fict - 1; i++){
			double dx = x[i + 1] - x[i];
			U_tilde[i][j] = U[i][j] - dt/dx*(F[i + 1][j] - F[i][j]); 
		}
	}
	BoundCond(U_tilde);

	DataArray W_tilde(Central_num);
	ConvertUtoW(W_tilde, U_tilde);
	//BoundCond(W_tilde);

	DataArray F_tilde(Central_num);
	InitialZeros(F_tilde, 3);

	for (int i = 0; i < Central_num; i++) {
		F_tilde[i][0] = W_tilde[i][0] * W_tilde[i][1]; // rho*u

		F_tilde[i][1] = W_tilde[i][0] * std::pow(W_tilde[i][1], 2) + W_tilde[i][2]; // rho*u^2 + P
		
		double E = 0.5 * std::pow(W_tilde[i][1], 2) + W_tilde[i][2] / (W_tilde[i][0] * (gamm - 1.0));
		F_tilde[i][2] = W_tilde[i][1] * (W_tilde[i][0] * E + W_tilde[i][2]); // u(rho*E + P)	
	}

	DataArray U_new(Central_num);
	InitialZeros(U_new, 3);

	for (int j = 0; j < 3; j++) {
		for (int i = fict; i < N + fict - 1; i++){
			double dx = x[i + 1] - x[i];
			U_new[i][j] = 0.5*(U[i][j] + U_tilde[i][j] - dt/dx*(F_tilde[i][j] - F_tilde[i - 1][j])); 
		}
	}
	BoundCond(U_new);
	ConvertUtoW(W_new, U_new);

	for (int i = fict; i < N + fict - 1; i++) {
		W[i] = W_new[i];
	}

	

}

void DOperator(
	std::vector<std::vector<double>>& DU,
	std::vector<std::vector<double>> U,
	double Q) {

	for (int j = 0; j < 3; j++) {
		for (int i = fict; i < N + fict - 1; i++) {
			DU[i][j] = Q * ((U[i + 1][j] - U[i][j]) - (U[i][j] - U[i - 1][j]));
		}
	}
}

void NOperator(
	std::vector<std::vector<double>>& NU, 
	std::vector<std::vector<double>> U,
	double Q) {

	std::vector<std::vector<double>> phi(N + 2 * fict);
	InitialZeros(phi, 3);
	std::vector<std::vector<double>> phi_c(N + 2 * fict);
	InitialZeros(phi_c, 3);
	std::vector<std::vector<double>> DU(N + 2 * fict - 1);
	InitialZeros(DU, 3);
	int s;

	DOperator(DU, U, Q);

	for (int j = 0; j < 3; j++) {
		for (int i = fict; i < N + fict; i++) {
			phi[i][j] = Q * (U[i][j] - U[i - 1][j]);
		}

		for (int i = fict; i < N + fict; i++) {
			if (phi[i][j] > 0) {
				s = 1;
			}
			else if (phi[i][j] < 0) {
				s = -1;
			}
			else {
				s = 0;
			}

			phi_c[i][j] = 
				s * std::max(0.0, 
					std::min(s * ((U[i - 1][j] + DU[i - 1][j]) - (U[i - 2][j] + DU[i - 2][j])), 
						std::min(std::abs(phi[i][j]),
							s * ((U[i + 1][j] + DU[i + 1][j]) - (U[i][j] + DU[i][j])))));
		}

		for (int i = fict; i < N + fict - 1; i++) {
			NU[i][j] = -(phi_c[i + 1][j] - phi_c[i][j]);
		}
	}	
}	



void UpdateArrays(
	std::vector<std::vector<double>>& W,
	std::vector<std::vector<double>> W_new,
	std::vector<std::vector<double>> W_b,
	std::vector<std::vector<double>> F,
	std::string method,
	std::string high_order_method,
	std::string solver,
	std::string TVD_solver,
	bool Viscous_flag, double mu0,
	std::string time_method,
	std::vector<double> x,double dt) {

	std::vector<std::vector<double>> W_L(N + 2 * fict);
	std::vector<std::vector<double>> W_R(N + 2 * fict);

	

	if (method == "MacCORMACK") {
		MacCORMACK(W, W_new, x, Viscous_flag, mu0, dt, solver);
		return;
	}
	
	for (int i = fict; i < N + fict - 1; i++) {
                double dx = x[i + 1] - x[i];
                double du = W[i + 1][1] - W[i][1];

                if (du / dx < 0.0) {
                        W[i][2] = W[i][2] + mu0 * W[i][0] * std::pow(du, 2);
                } else {
                        W[i][2] = W[i][2];
                }
	}

	if (time_method == "RK3") {
		RK3(W_new, W, method, solver, fict, N + fict - 1, x, dt);
	} else {
		Euler(W_new, W, method, high_order_method, solver, TVD_solver, fict, N + fict - 1, x, dt);
	}	

	for (int i = fict; i < N + fict - 1; i++) {
		W[i] = W_new[i];
	}
}
