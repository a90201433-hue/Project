#include <vector>
#include <cmath>
#include <iostream>
#include "RiemannSolver.h"
#include "Types.h"

extern int Nx, Ny, fict;
extern double gamm;

void ConvertWtoU(const Field W, Field& U, int dir) {
    int nx = W.size();
    U.resize(nx);
    
    for (int i = 0; i < nx; i++) {
        int ny = W[i].size();
        U[i].resize(ny);
        for (int j = 0; j < ny; j++) {
            double rho = W[i][j][0];
            double u   = W[i][j][1];
            double v   = W[i][j][2];
            double P   = W[i][j][3];

            double E = P / (rho * (gamm - 1.0)) + 0.5 * (u * u + v * v);

            U[i][j][0] = rho;
            U[i][j][1] = rho * u;
            U[i][j][2] = rho * v;
            U[i][j][3] = rho * E;
        }
    }
}

void ConvertUtoW(Field& W, const Field U, int dir) {
    int nx = U.size();
    W.resize(nx);

    for (int i = 0; i < nx; i++) {
        int ny = U[i].size();
        W[i].resize(ny);
        
        for (int j = 0; j < ny; j++) {
            double rho = std::max(1e-6, U[i][j][0]);
            double u   = U[i][j][1] / rho;
            double v   = U[i][j][2] / rho;
            double E   = U[i][j][3] / rho;
            
            double P = std::max(1e-6, (gamm - 1.0) * rho * (E - 0.5 * (u*u + v*v)));
            
            W[i][j][0] = rho;
            W[i][j][1] = u;
            W[i][j][2] = v;
            W[i][j][3] = P;
        }
    }
}

void ConvertWtoU(State W, State& U) {
    
	double rho = W[0];
	double u   = W[1];
	double v   = W[2];
	double P   = W[NEQ - 1];

	double E = P / (rho * (gamm - 1.0)) 
				+ 0.5 * (u * u + v * v);

	U[0] = rho;
	U[1] = rho * u;
	U[2] = rho * v;
	U[3] = rho * E;

	return;	
}
