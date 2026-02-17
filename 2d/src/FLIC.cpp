#include <cmath>
#include <vector>
#include "FLIC.h"
#include "Init.h"
#include "Types.h"

extern double gamm, Lx, Ly, C1, C2;
extern int Nx, Ny, fict;

void FLIC_L(Field W,
            Field& W_tilde,
            std::vector<double> x,
            std::vector<double> y,
            double dt,
            int dir) {

    if (dir == 0) {
        for (int i = fict; i < Nx + fict - 1; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                double dx = x[i + 1] - x[i];

                double rho = W[i][j][0];
                double u = W[i][j][1];
                double v = W[i][j][2];
                double P = W[i][j][NEQ - 1];

                double P_ip = 0.5 * (W[i][j][NEQ - 1] + W[i + 1][j][NEQ - 1]);
                double P_im = 0.5 * (W[i - 1][j][NEQ - 1] + W[i][j][NEQ - 1]);

                W_tilde[i][j][0] = rho;

                double u_tilde = u - dt / rho * (P_ip - P_im) / dx;
                W_tilde[i][j][1] = u_tilde;

                W_tilde[i][j][2] = v;

                double E = P / (gamm - 1.0) + 0.5 * rho * (u * u + v * v);
                double u_ip = 0.5 * (W[i][j][1] + W[i + 1][j][1]);
                double u_im = 0.5 * (W[i - 1][j][1] + W[i][j][1]);
                double E_tilde = E - dt * (P_ip * u_ip - P_im * u_im) / dx;
                double kinetic_tilde = 0.5 * rho * (u_tilde * u_tilde + v * v);
                double P_tilde = (gamm - 1.0) * (E_tilde - kinetic_tilde);
        
                W_tilde[i][j][NEQ - 1] = std::max(1e-8, P_tilde);
            }
        }
    }

    if (dir == 1) {
        for (int i = fict; i < Nx + fict; i++) {
            for (int j = fict; j < Ny + fict - 1; j++) {
                double dy = y[j + 1] - y[j];

                double rho = W[i][j][0];
                double u = W[i][j][1];
                double v = W[i][j][2];
                double P = W[i][j][NEQ - 1];

                double P_jp = 0.5 * (W[i][j][NEQ - 1] + W[i][j + 1][NEQ - 1]);
                double P_jm = 0.5 * (W[i][j - 1][NEQ - 1] + W[i][j][NEQ - 1]);

                W_tilde[i][j][0] = rho;

                W_tilde[i][j][1] = u;

                double v_tilde = v - dt / rho * (P_jp - P_jm) / dy;
                W_tilde[i][j][2] = v_tilde;

                double E = P / (gamm - 1.0) + 0.5 * rho * (u * u + v * v);
                double v_jp = 0.5 * (W[i][j][2] + W[i][j + 1][2]);
                double v_jm = 0.5 * (W[i][j - 1][2] + W[i][j][2]);
                double E_tilde = E - dt  * (P_jp * v_jp - P_jm * v_jm) / dy;
                double kinetic_tilde = 0.5 * rho * (u * u + v_tilde * v_tilde);
                double P_tilde = (gamm - 1.0) * (E_tilde - kinetic_tilde);

                W_tilde[i][j][NEQ - 1] = std::max(1e-8, P_tilde);
            }
        }
    }
}

void FLIC_E(Field W_tilde,
            Field& W_new,
            std::vector<double> x,
            std::vector<double> y,
            double dt,
            int dir) {

    if (dir == 0) {  
        for (int i = fict; i < Nx + fict - 1; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                double dx = x[i + 1] - x[i];

                // upwind на границе i+1/2 
                double u_face = 0.5 * (W_tilde[i][j][1] + W_tilde[i+1][j][1]);

                int up = (u_face > 0.0) ? i : i+1;

                double rho_flux  = W_tilde[up][j][0] * W_tilde[up][j][1];
                double momx_flux = rho_flux * W_tilde[up][j][1];
                double momy_flux = rho_flux * W_tilde[up][j][2];

                double E = W_tilde[up][j][NEQ - 1] / (gamm - 1.0) + 
                            0.5 * W_tilde[up][j][0] * 
                            (W_tilde[up][j][1]*W_tilde[up][j][1] + W_tilde[up][j][2] * W_tilde[up][j][2]);

                double E_flux = W_tilde[up][j][1] * E;

                // аналогично для i-1/2 
                double u_face_m = 0.5 * (W_tilde[i - 1][j][1] + W_tilde[i][j][1]);

                int up_m = (u_face_m > 0.0) ? i-1 : i;

                double rho_flux_m  = W_tilde[up_m][j][0] * W_tilde[up_m][j][1];
                double momx_flux_m = rho_flux_m * W_tilde[up_m][j][1];
                double momy_flux_m = rho_flux_m * W_tilde[up_m][j][2];

                double E_m = W_tilde[up_m][j][NEQ - 1] / (gamm - 1.0) + 0.5 * W_tilde[up_m][j][0] * (W_tilde[up_m][j][1]*W_tilde[up_m][j][1] + W_tilde[up_m][j][2]*W_tilde[up_m][j][2]);

                double E_flux_m = W_tilde[up_m][j][1] * E_m;

                // обновление 
                double rho_new = W_tilde[i][j][0] - dt/dx * (rho_flux - rho_flux_m);

                double momx_new = W_tilde[i][j][0] * W_tilde[i][j][1] - dt/dx * (momx_flux - momx_flux_m);

                double momy_new = W_tilde[i][j][0] * W_tilde[i][j][2] - dt/dx * (momy_flux - momy_flux_m);

                double E_new = (W_tilde[i][j][NEQ - 1] / (gamm - 1.0) + 0.5 * W_tilde[i][j][0] * (W_tilde[i][j][1]*W_tilde[i][j][1] + W_tilde[i][j][2]*W_tilde[i][j][2])) - dt/dx * (E_flux - E_flux_m);

                rho_new = std::max(1e-8, rho_new);

                double u_new = momx_new / rho_new;
                double v_new = momy_new / rho_new;

                double kinetic = 0.5 * rho_new * (u_new * u_new + v_new * v_new);
                double P_new = (gamm - 1.0) * (E_new - kinetic);

                W_new[i][j][0] = rho_new;
                W_new[i][j][1] = u_new;
                W_new[i][j][2] = v_new;
                W_new[i][j][3] = std::max(1e-8, P_new);
            }
        }
    }

    if (dir == 1) {   
        for (int i = fict; i < Nx + fict; i++) {
            for (int j = fict; j < Ny + fict - 1; j++) {
                double dy = y[j + 1] - y[j];

                double v_face = 0.5 * (W_tilde[i][j][2] + W_tilde[i][j+1][2]);
                int up = (v_face > 0.0) ? j : j+1;

                double rho_flux  = W_tilde[i][up][0] * W_tilde[i][up][2];
                double momx_flux = rho_flux * W_tilde[i][up][1];
                double momy_flux = rho_flux * W_tilde[i][up][2];

                double E = W_tilde[i][up][3] / (gamm - 1.0) + 0.5 * W_tilde[i][up][0] * (W_tilde[i][up][1]*W_tilde[i][up][1] + W_tilde[i][up][2]*W_tilde[i][up][2]);

                double E_flux = W_tilde[i][up][2] * E;

                double v_face_m = 0.5 * (W_tilde[i][j-1][2] + W_tilde[i][j][2]);
                int up_m = (v_face_m > 0.0) ? j-1 : j;

                double rho_flux_m  = W_tilde[i][up_m][0] * W_tilde[i][up_m][2];
                double momx_flux_m = rho_flux_m * W_tilde[i][up_m][1];
                double momy_flux_m = rho_flux_m * W_tilde[i][up_m][2];

                double E_m = W_tilde[i][up_m][3] / (gamm - 1.0) + 0.5 * W_tilde[i][up_m][0] * (W_tilde[i][up_m][1]*W_tilde[i][up_m][1] + W_tilde[i][up_m][2]*W_tilde[i][up_m][2]);

                double E_flux_m = W_tilde[i][up_m][2] * E_m;

                double rho_new = W_tilde[i][j][0] - dt/dy * (rho_flux - rho_flux_m);

                double momx_new = W_tilde[i][j][0] * W_tilde[i][j][1] - dt/dy * (momx_flux - momx_flux_m);

                double momy_new = W_tilde[i][j][0] * W_tilde[i][j][2] - dt/dy * (momy_flux - momy_flux_m);

                double E_new = (W_tilde[i][j][3] / (gamm - 1.0) + 0.5 * W_tilde[i][j][0] * (W_tilde[i][j][1]*W_tilde[i][j][1] + W_tilde[i][j][2]*W_tilde[i][j][2])) - dt/dy * (E_flux - E_flux_m);

                rho_new = std::max(1e-8, rho_new);

                double u_new = momx_new / rho_new;
                double v_new = momy_new / rho_new;

                double kinetic = 0.5 * rho_new * (u_new * u_new + v_new * v_new);
                double P_new = (gamm - 1.0) * (E_new - kinetic);

                W_new[i][j][0] = rho_new;
                W_new[i][j][1] = u_new;
                W_new[i][j][2] = v_new;
                W_new[i][j][3] = std::max(1e-8, P_new);
            }
        }
    }
}




