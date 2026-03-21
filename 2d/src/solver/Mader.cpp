#include <cmath>
#include <vector>
#include <iostream>
#include "Mader.h"
#include "FileProcessing.h"
#include "BoundCond.h"
#include "Init.h"
#include "Types.h"

extern double gamm, Lx, Ly, T_init, R_gas, M, P_min, E_act, Z_freq, VISC, MINWT, GASW, MINGRHO;
extern int Nx, Ny, fict;

void EOS(Field rho, Field I, Field& P, Field& T) {
    // Тут проблема, в методе Мейдера на вход поступает внутренняя энергия, а не давление
    // Если мы будем считать из W внутреннюю энергию, а из нее давление
    // То получится масло масялное, давление будет вычисляться само через себя
    // Хз пока, что с этим делать, поэтому я пока оно само через себя
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            P[i][j][0] = (((gamm - 1.0) * rho[i][j][0] * I[i][j][0]) < P_min) ? P_min : ((gamm - 1.0) * rho[i][j][0] * I[i][j][0]);
            T[i][j][0] = P[i][j][0] / (rho[i][j][0] * R_gas / M);
        }
    }
}

void Arrenius(Field& mass_fraction, Field T, double dt, int Nx_tot, int Ny_tot) {
    // Здесь я пока полностью заигнорила проверку условия на ложную детонацию, которая зависит от количества
    // пройденных циклов по времени
    Field mass_fraction_new(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            mass_fraction_new[i][j][0] = mass_fraction[i][j][0] - dt * Z_freq * mass_fraction[i][j][0] * exp(-E_act / (R_gas * T[i][j][0]));
            mass_fraction_new[i][j][0] = (mass_fraction_new[i][j][0] < GASW && T[i][j][0] < MINWT) ? 0.0 : mass_fraction_new[i][j][0];
        }
    }
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            mass_fraction[i][j][0] = mass_fraction_new[i][j][0];
        }
    }
}

void Viscosity(Field& q1, Field& q2, Field& q3, Field& q4, Field rho, Field u, Field v) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            q1[i][j][0] = (v[i][j][0] >= v[i][j + 1][0]) ? VISC * rho[i][j][0] * (v[i][j][0] - v[i][j + 1][0]) : 0.0;
            q2[i][j][0] = (u[i][j][0] >= u[i + 1][j][0]) ? VISC * rho[i][j][0] * (u[i][j][0] - u[i + 1][j][0]) : 0.0;
        }
    }
    for (int i = fict; i < Nx + fict - 1; i++)
    for (int j = fict; j < Ny + fict - 1; j++) {
        q3[i][j][0] = (j > fict) ? q1[i][j-1][0] : 0.0;
        q4[i][j][0] = (i > fict) ? q2[i-1][j][0] : 0.0;
    }
}

void VelocityTilde(Field& u_tilde, Field& v_tilde, Field P, Field rho, const std::vector<double>& x, const std::vector<double>& y, Field q1, Field q2, Field q3, Field q4, Field u, Field v, double dt) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double P1 = P[i][j - 1][0];
            double P2 = P[i - 1][j][0];
            double P3 = P[i][j + 1][0];
            double P4 = P[i + 1][j][0];

            double dx = x[i] - x[i - 1];
            double dy = y[j] - y[j - 1];

            v_tilde[i][j][0] = v[i][j][0] - dt / (rho[i][j][0] * 2.0 * dy) * ((P3 - P1) + (q3[i][j][0] - q1[i][j][0]));
            u_tilde[i][j][0] = u[i][j][0] - dt / (rho[i][j][0] * 2.0 * dx) * ((P4 - P2) + (q4[i][j][0] - q2[i][j][0]));
        }
    }
}

void ZIPEnergy(Field I, Field& I_tilde, Field P, double dt, Field rho, Field u, Field u_tilde, Field v, Field v_tilde, Field q1, Field q2, Field q3, Field q4, const std::vector<double>& x, const std::vector<double>& y) {
    // Здесь только уравнение сосотояния идеального газа для I и P!
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double U1 = u[i - 1][j][0] + u_tilde[i - 1][j][0];
            double U2 = u[i + 1][j][0] + u_tilde[i + 1][j][0];
            double V1 = v[i][j - 1][0] + v_tilde[i][j - 1][0];
            double V2 = v[i][j + 1][0] + v_tilde[i][j + 1][0];
            double T3 = u[i][j][0] + u_tilde[i][j][0];
            double T1 = v[i][j][0] + v_tilde[i][j][0];
            double dx = x[i] - x[i - 1];
            double dy = y[j] - y[j - 1];

            I_tilde[i][j][0] = I[i][j][0] - dt / (4 * rho[i][j][0]) * ((P[i][j][0] / dx) * (U2 - U1) +
                                                                        (q4[i][j][0] / dx) * (U2 - T3) +
                                                                        (q2[i][j][0] / dx) * (T3 - U1) +
                                                                        (P[i][j][0] / dy) * (V2 - V1) +
                                                                        (q3[i][j][0] / dy) * (V2 - T1) +
                                                                        (q1[i][j][0] / dy) * (T1 - V1));
        }
    }
}

void ChangingFluxes(Field& alpha, Field& beta, Field& DM, Field& DE, Field& DW, Field& DPU, Field& DPV, Field u_tilde, Field v_tilde, const std::vector<double>& x, const std::vector<double>& y, double dt, Field rho, Field I_tilde, Field u, Field v, Field mass_fraction) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double dx = x[i] - x[i - 1];
            double dy = y[j] - y[j - 1];

            alpha[i][j][0] = (0.5 * (u_tilde[i - 1][j][0] + u_tilde[i][j][0]) * dt / dx) / (1 + (u_tilde[i - 1][j][0] - u_tilde[i][j][0]) * dt / dx);
            beta[i][j][0] = (0.5 * (v_tilde[i][j - 1][0] + v_tilde[i][j][0]) * dt / dy) / (1 + (v_tilde[i][j - 1][0] - v_tilde[i][j][0]) * dt / dy);
 
            if (alpha[i][j][0] >= 0) {
                double DMASS = rho[i - 1][j][0] * std::abs(alpha[i][j][0]);
 
                DM[i][j][0] += DMASS;
                DM[i - 1][j][0] -= DMASS;

                DE[i][j][0] += DMASS * (I_tilde[i - 1][j][0] + 0.5 * (u_tilde[i - 1][j][0] * u_tilde[i - 1][j][0] + v_tilde[i - 1][j][0] * v_tilde[i - 1][j][0]));
                DE[i - 1][j][0] -= DMASS * (I_tilde[i - 1][j][0] + 0.5 * (u_tilde[i - 1][j][0] * u_tilde[i - 1][j][0] + v_tilde[i - 1][j][0] * v_tilde[i - 1][j][0]));
            
                DPU[i][j][0] += DMASS * u_tilde[i - 1][j][0];
                DPU[i - 1][j][0] -= DMASS * u_tilde[i - 1][j][0];

                DPV[i][j][0] += DMASS * v_tilde[i - 1][j][0];
                DPV[i - 1][j][0] -= DMASS * v_tilde[i - 1][j][0];
            
                DW[i][j][0] += DMASS * mass_fraction[i - 1][j][0];
                DW[i - 1][j][0] -= DMASS * mass_fraction[i - 1][j][0];
            }
            if (alpha[i][j][0] < 0) {
                double DMASS = rho[i][j][0] * std::abs(alpha[i][j][0]);

                DM[i][j][0] -= DMASS;
                DM[i - 1][j][0] += DMASS;

                DE[i][j][0] -= DMASS * (I_tilde[i][j][0] + 0.5 * (u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]));
                DE[i - 1][j][0] += DMASS * (I_tilde[i][j][0] + 0.5 * (u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]));
            
                DPU[i][j][0] -= DMASS * u_tilde[i][j][0];
                DPU[i - 1][j][0] += DMASS * u_tilde[i][j][0];

                DPV[i][j][0] -= DMASS * v_tilde[i][j][0];
                DPV[i - 1][j][0] += DMASS * v_tilde[i][j][0];
            
                DW[i][j][0] -= DMASS * mass_fraction[i][j][0];
                DW[i - 1][j][0] += DMASS * mass_fraction[i][j][0];
            }
            if (beta[i][j][0] >= 0) {
                double DMASS = rho[i][j - 1][0] * std::abs(beta[i][j][0]);

                DM[i][j][0] += DMASS;
                DM[i][j - 1][0] -= DMASS;
            
                DE[i][j][0] += DMASS * (I_tilde[i][j - 1][0] + 0.5 * (u_tilde[i][j - 1][0] * u_tilde[i][j - 1][0] + v_tilde[i][j - 1][0] * v_tilde[i][j - 1][0]));
                DE[i][j - 1][0] -= DMASS * (I_tilde[i][j - 1][0] + 0.5 * (u_tilde[i][j - 1][0] * u_tilde[i][j - 1][0] + v_tilde[i][j - 1][0] * v_tilde[i][j - 1][0]));
            
                DPU[i][j][0] += DMASS * u_tilde[i][j - 1][0];
                DPU[i][j - 1][0] -= DMASS * u_tilde[i][j - 1][0];

                DPV[i][j][0] += DMASS * v_tilde[i][j - 1][0];
                DPV[i][j - 1][0] -= DMASS * v_tilde[i][j - 1][0];
            
                DW[i][j][0] += DMASS * mass_fraction[i][j - 1][0];
                DW[i][j - 1][0] -= DMASS * mass_fraction[i][j - 1][0];
            }
            if (beta[i][j][0] < 0) {
                double DMASS = rho[i][j][0] * std::abs(beta[i][j][0]);

                DM[i][j][0] -= DMASS;
                DM[i][j - 1][0] += DMASS;
            
                DE[i][j][0] -= DMASS * (I_tilde[i][j][0] + 0.5 * (u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]));
                DE[i][j - 1][0] += DMASS * (I_tilde[i][j][0] + 0.5 * (u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]));
            
                DPU[i][j][0] -= DMASS * u_tilde[i][j][0];
                DPU[i][j - 1][0] += DMASS * u_tilde[i][j][0];

                DPV[i][j][0] -= DMASS * v_tilde[i][j][0];
                DPV[i][j - 1][0] += DMASS * v_tilde[i][j][0];
            
                DW[i][j][0] -= DMASS * mass_fraction[i][j][0];
                DW[i][j - 1][0] += DMASS * mass_fraction[i][j][0];
            }
        }
    }
    for (int j = fict; j < Ny + fict - 1; j++) {
        int i = Nx + fict - 2; 
        double dx = x[i] - x[i-1];
        
        double alpha_right = u_tilde[i][j][0] * dt / dx;
        
        if (alpha_right > 0) { 
            double DMASS = rho[i][j][0] * alpha_right;
            double E_d = I_tilde[i][j][0] + 0.5*(u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]);
            DM[i][j][0]  -= DMASS;
            DE[i][j][0]  -= E_d * DMASS;
            DPU[i][j][0] -= u_tilde[i][j][0] * DMASS;
            DPV[i][j][0] -= v_tilde[i][j][0] * DMASS;
            DW[i][j][0]  -= mass_fraction[i][j][0] * DMASS;
        }
    }
    for (int i = fict; i < Nx + fict - 1; i++) {
        int j = Ny + fict - 2;
        double dy = y[j] - y[j-1];
        
        double beta_top = v_tilde[i][j][0] * dt / dy;
        
        if (beta_top > 0) {
            double DMASS = rho[i][j][0] * beta_top;
            double E_d = I_tilde[i][j][0] + 0.5*(u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]);
            DM[i][j][0]  -= DMASS;
            DE[i][j][0]  -= E_d * DMASS;
            DPU[i][j][0] -= u_tilde[i][j][0] * DMASS;
            DPV[i][j][0] -= v_tilde[i][j][0] * DMASS;
            DW[i][j][0]  -= mass_fraction[i][j][0] * DMASS;
        }
    }
}

void Repartition(Field& rho, Field& u, Field& v, Field& I, Field& mass_fraction, Field I_tilde, Field u_tilde, Field v_tilde, Field DM, Field DE, Field DW, Field DPU, Field DPV) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double rho_new = rho[i][j][0] + DM[i][j][0];

            if (rho_new <= 1e-12) {
                rho[i][j][0] = 0.0;
                u[i][j][0] = 0.0;
                v[i][j][0] = 0.0;
                I[i][j][0] = 0.0;
                mass_fraction[i][j][0] = 0.0;
                continue;
            }

            double u_new = (rho[i][j][0] * u_tilde[i][j][0] + DPU[i][j][0]) / rho_new;
            double v_new = (rho[i][j][0] * v_tilde[i][j][0] + DPV[i][j][0]) / rho_new;
            double E = I_tilde[i][j][0] + 0.5 * (u_tilde[i][j][0] * u_tilde[i][j][0] + v_tilde[i][j][0] * v_tilde[i][j][0]);
            double I_new = (rho[i][j][0] * E + DE[i][j][0]) / rho_new - 0.5 * (u_new * u_new + v_new * v_new);
            double mass_fraction_new = (rho[i][j][0] * mass_fraction[i][j][0] + DW[i][j][0]) / rho_new;

            rho[i][j][0] = rho_new;
            u[i][j][0] = u_new;
            v[i][j][0] = v_new;
            I[i][j][0] = I_new;
            mass_fraction[i][j][0] = mass_fraction_new;
        }
    }
}


void Mader(Field& W_new, const Field& W, const std::vector<double>& x, const std::vector<double>& y, double dt) {

    int Nx_tot = Nx + 2*fict - 1;
    int Ny_tot = Ny + 2*fict - 1;

    Field mass_fraction(Nx_tot, std::vector<State>(Ny_tot, {1.0, 0.0, 0.0, 0.0}));
    Field rho(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            rho[i][j] = {};
            rho[i][j][0] = W[i][j][0];
        }
    }
    Field P(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            P[i][j] = {};
            P[i][j][0] = W[i][j][NEQ - 1];
        }
    }
    Field I(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            I[i][j] = {};
            I[i][j][0] = W[i][j][NEQ-1] / ((gamm - 1.0) * W[i][j][0]);
        }
    }
    Field T(Nx_tot, std::vector<State>(Ny_tot, {T_init, 0.0, 0.0, 0.0}));
    Field q1(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field q2(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field q3(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field q4(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field u(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            u[i][j] = {};
            u[i][j][0] = W[i][j][1];
        }
    }
    Field v(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            v[i][j] = {};
            v[i][j][0] = W[i][j][2];
        }
    }
    Field u_tilde(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            u_tilde[i][j] = {};
            u_tilde[i][j][0] = W[i][j][1];
        }
    }
    Field v_tilde(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            v_tilde[i][j] = {};
            v_tilde[i][j][0] = W[i][j][2];
        }
    }
    Field I_tilde(Nx_tot, std::vector<State>(Ny_tot));
    for (int i = 0; i < W.size(); i++) {
        for (int j = 0; j < W[i].size(); j++) {
            I_tilde[i][j] = {};
        }
    }
    Field alpha(Nx_tot, std::vector<State>(Ny_tot));
    Field beta(Nx_tot, std::vector<State>(Ny_tot));

    Field DM(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DE(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DW(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DPU(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DPV(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));

    EOS(rho, I, P, T);
    BoundCond(P);
    BoundCond(T);

    Arrenius(mass_fraction, T, dt, Nx_tot, Ny_tot);
    BoundCond(mass_fraction);

    Viscosity(q1, q2, q3, q4, rho, u, v);
    BoundCond(q1);
    BoundCond(q2);
    BoundCond(q3);
    BoundCond(q4);

    VelocityTilde(u_tilde, v_tilde, P, rho, x, y, q1, q2, q3, q4, u, v, dt);
    BoundCond(u_tilde);
    BoundCond(v_tilde);

    ZIPEnergy(I, I_tilde, P, dt, rho, u, u_tilde, v, v_tilde, q1, q2, q3, q4, x, y);
    BoundCond(I_tilde);

    ChangingFluxes(alpha, beta, DM, DE, DW, DPU, DPV, u_tilde, v_tilde, x, y, dt, rho, I_tilde, u, v, mass_fraction);

    Repartition(rho, u, v, I, mass_fraction, I_tilde, u_tilde, v_tilde, DM, DE, DW, DPU, DPV);
    BoundCond(rho);
    BoundCond(u);
    BoundCond(v);
    BoundCond(I);
    BoundCond(mass_fraction);
    
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {

            W_new[i][j][0] = rho[i][j][0];
            W_new[i][j][1] = u[i][j][0];
            W_new[i][j][2] = v[i][j][0];

            double kinetic = 0.5 * rho[i][j][0] * (u[i][j][0] * u[i][j][0] + v[i][j][0] * v[i][j][0]);
            double E = rho[i][j][0] * I[i][j][0] + kinetic;

            double P_new = (gamm - 1.0) * (E - kinetic);

            W_new[i][j][NEQ - 1] = std::max(P_min, P_new);
        }
    }
}