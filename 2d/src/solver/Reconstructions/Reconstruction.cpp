#include <vector>
#include "Limiters.h"
#include "Types.h"

//extern double gamm, Lx, Ly, C1, C2;
extern int Nx, Ny, fict;
extern std::string method, rec_limiter;

void ReconstructGodunov(const Field& W,
                        Field& W_L,
                        Field& W_R,
                        int dir) {
    if (dir == 0) {
        for (int j = fict; j < Ny - 1 + fict; j++) {
            for (int i = fict; i < Nx + fict; i++) {
                W_L[i][j] = W[i - 1][j]; // левая ячейка
                W_R[i][j] = W[i][j];     // правая ячейка
            }
        }
    }

    else if (dir == 1) {
        for (int i = fict; i < Nx - 1 + fict; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                W_L[i][j] = W[i][j - 1]; // нижняя ячейка
                W_R[i][j] = W[i][j];     // верхняя ячейка
            }
        }
    }
}


void ReconstructKolgan(const Field& W,
                       Field& W_L,
                       Field& W_R,
                       int dir) {
    if (dir == 0) {
        for (int j = fict; j < Ny - 1 + fict; j++) {
            for (int i = fict; i < Nx + fict; i++) {
                State dWm, dWp;

                for (int k = 0; k < NEQ; k++) {
                    dWm[k] = W[i-1][j][k] - W[i-2][j][k];
                    dWp[k] = W[i][j][k]   - W[i-1][j][k];
                }

                State slope;
                if (rec_limiter == "minmod")
                    slope = Minmod(dWm, dWp);
                else if (rec_limiter == "superbee")
                    slope = Superbee(dWm, dWp);
                else if (rec_limiter == "vanleer")
                    slope = Vanleer(dWm, dWp);

                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W[i-1][j][k] + 0.5 * slope[k];

                    W_R[i][j][k] =
                        W[i][j][k]   - 0.5 * slope[k];
                }
            }
        }
    }

    else if (dir == 1) {
        for (int i = fict; i < Nx - 1 + fict; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                State dWm, dWp;

                for (int k = 0; k < NEQ; k++) {
                    dWm[k] = W[i][j-1][k] - W[i][j-2][k];
                    dWp[k] = W[i][j][k]   - W[i][j-1][k];
                }

                State slope;
                if (rec_limiter == "minmod")
                    slope = Minmod(dWm, dWp);
                else if (rec_limiter == "superbee")
                    slope = Superbee(dWm, dWp);
                else if (rec_limiter == "vanleer")
                    slope = Vanleer(dWm, dWp);

                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W[i][j-1][k] + 0.5 * slope[k];

                    W_R[i][j][k] =
                        W[i][j][k]   - 0.5 * slope[k];
                }
            }
        }
    }
}

void Reconstruct(const Field& W,
                 Field& W_L,
                 Field& W_R,
                 int dir) {
	
	if (method == "Godunov")
        ReconstructGodunov(W, W_L, W_R, dir);

    else if (method == "Kolgan")
        ReconstructKolgan(W, W_L, W_R, dir);

	return;
}