#include <vector>
#include <cmath>
#include <iostream>

extern int N, fict;
extern double gamm;

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

void ConvertUtoW(std::vector<std::vector<double>>& W,
		 std::vector<std::vector<double>> U){

	int K = U.size();
	for (int i = 0; i < K; i ++){
		//W[i] = {0, 0, 0};
		W[i].resize(3);

		W[i][0] = U[i][0]; // rho

		W[i][1] = U[i][1]/U[i][0]; // u

		W[i][2] = (gamm - 1)*(U[i][2] - 0.5*W[i][0]*std::pow(W[i][1], 2)); // P
	}
	
	for (int i = 0; i < K; i ++){
		if (W[i][0] < 1e-10){
			W[i][0] = 1e-10;
		}

		if (W[i][2] < 1e-6){
			W[i][2] = 1e-10;
		}
	}
}
