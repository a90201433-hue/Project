#include <vector>
#include <cmath>
#include <iostream>

extern int N, fict;
extern double gamm;

using DataArray = std::vector<std::vector<double>>;

void ConvertWtoU(std::vector<std::vector<double>> W,
		 std::vector<std::vector<double>>& U){

	int K = W.size();
	double E = 0.0;
	for (int i = 0; i < K; i++){
		//U[i] = {0, 0, 0};
		U[i].resize(3);

		U[i][0] = W[i][0]; // rho

		U[i][1] = W[i][0] * W[i][1]; // rho*u

		E = 0.5*std::pow(W[i][1], 2) + W[i][2]/(W[i][0]*(gamm - 1));
		U[i][2] = W[i][0]*E ; // rho*E
	}
}

void ConvertWtoU(std::vector<double> W,
		 std::vector<double>& U) {

	double E = 0.0;
	U = {0, 0, 0};

	U[0] = W[0]; // rho
	U[1] = W[0] * W[1]; // rho*u
	E = 0.5*std::pow(W[1], 2) + W[2]/(W[0]*(gamm - 1));
	U[2] = W[0]*E ;
	
	return;
}

void ConvertUtoW(std::vector<std::vector<double>>& W,
		 std::vector<std::vector<double>> U){

	int K = U.size();
	for (int i = 0; i < K; i ++){
		//W[i] = {0, 0, 0};
		W[i].resize(3);

		W[i][0] = std::max(1e-6, U[i][0]); // rho

		W[i][1] = U[i][1]/U[i][0]; // u

		W[i][2] = std::max(1e-6, (gamm - 1)*(U[i][2] - 0.5*W[i][0]*std::pow(W[i][1], 2))); // P
	}
	
	/*for (int i = 0; i < K; i ++){
		if (W[i][0] < 1e-6){
			W[i][0] = 1e-6;
		}

		if (W[i][2] < 1e-6){
			W[i][2] = 1e-6;
		}
	}*/
}

void ConvertUtoW(std::vector<double>& W,
		 std::vector<double> U) {

	W = {0, 0, 0};

	W[0] = std::max(1e-6, U[0]);
	W[1] = U[1]/U[0];
	W[2] = std::max(1e-6, (gamm - 1)*(U[2] - 0.5*W[0]*std::pow(W[1], 2)));
/*
	if (W[0] < 1e-6) {
		W[0] = 1e-6;
	}

	if (W[2] < 1e-6) {
		W[2] = 1e-6;
	}
	*/
	return;
}

void FindFfromW(std::vector<double>& F, 
		std::vector<double> W) {
	
	F = {0, 0, 0};
		
	F[0] = W[0] * W[1]; // rho*u

	F[1] = W[0] * std::pow(W[1], 2) + W[2]; // rho*u^2 + P
		
	double E = 0.5 * std::pow(W[1], 2) + W[2] / (W[0] * (gamm - 1.0));
	F[2] = W[1] * (W[0] * E + W[2]);
}
