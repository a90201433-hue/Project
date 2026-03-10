#include <vector>
#include <cmath>
#include <iostream>
#include "RiemannSolver.h"
#include "Types.h"

extern int Nx, Ny, fict;
extern double gamma_gas;


/* ===============================
   W -> U  (Field)
   =============================== */
void ConvertWtoU(const Field W, Field& U, int dir)
{
    int nx = W.size();
    U.resize(nx);

    for (int i = 0; i < nx; i++)
    {
        int ny = W[i].size();
        U[i].resize(ny);

        for (int j = 0; j < ny; j++)
        {
            double rho = W[i][j][0];
            double u   = W[i][j][1];
            double v   = W[i][j][2];
            double P   = W[i][j][3];

            double E =
                P/(rho*(gamma_gas-1.0))
                +0.5*(u*u+v*v);

            U[i][j][0]=rho;
            U[i][j][1]=rho*u;
            U[i][j][2]=rho*v;
            U[i][j][3]=rho*E;
        }
    }
}


/* ===============================
   W -> U  (State)
   =============================== */

void ConvertWtoU(const State W, State& U)
{
    double rho = W[0];
    double u   = W[1];
    double v   = W[2];
    double P   = W[3];

    double E =
        P / (rho * (gamma_gas - 1.0))
        + 0.5 * (u*u + v*v);

    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = rho * E;
}


/* ===============================
   U -> W  (Field)
   =============================== */

void ConvertUtoW(Field& W, const Field U, int dir)
{
    for (int i = fict; i < Nx + fict - 1; i++)
    for (int j = fict; j < Ny + fict - 1; j++)
    {
        double rho = std::max(1e-8,U[i][j][0]);

        double u = U[i][j][1]/rho;
        double v = U[i][j][2]/rho;

        double E = U[i][j][3]/rho;

        double P =
        std::max(
        1e-8,
        (gamma_gas-1.0)*
        rho*
        (E-0.5*(u*u+v*v))
        );

        W[i][j][0]=rho;
        W[i][j][1]=u;
        W[i][j][2]=v;
        W[i][j][3]=P;
    }
}


/* ===============================
   U -> W  (State)
   =============================== */

void ConvertUtoW(const State U, State& W)
{
    double rho = std::max(1e-8, U[0]);

    double u = U[1] / rho;
    double v = U[2] / rho;

    double E = U[3] / rho;

    double P =
        std::max(
            1e-8,
            (gamma_gas - 1.0) *
            rho *
            (E - 0.5*(u*u + v*v))
        );

    W[0] = rho;
    W[1] = u;
    W[2] = v;
    W[3] = P;
}