#include <vector>
#include <string>
#include <iostream>
#include "Types.h"

extern int Nx, Ny, fict;
extern std::string x_left_bound, x_right_bound,
				   y_up_bound, y_down_bound;

void BoundCond(Field& W) {

	size_t Nx_cells = Nx - 1;
    size_t Ny_cells = Ny - 1;

	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

	
    // ---- X-границы ----
    for (size_t j = fict; j < Ny_tot - fict; j++) {
        for (size_t g = 0; g < fict; g++) {

            // ---- Left ----
            if (x_left_bound == "free")
                W[g][j] = W[fict][j];

            else if (x_left_bound == "wall") {
                W[g][j] = W[fict][j];
                W[g][j][1] *= -1.0;  // меняем u
            }

            else if (x_left_bound == "periodic")
                W[g][j] = W[Nx_cells + g][j];

            // ---- right ----
            if (x_right_bound == "free")
                W[Nx_cells + fict + g][j] =
                    W[Nx_cells + fict - 1][j];

            else if (x_right_bound == "wall") {
                W[Nx_cells + fict + g][j] =
                    W[Nx_cells + fict - 1][j];
                W[Nx_cells + fict + g][j][1] *= -1.0;
            }

            else if (x_right_bound == "periodic")
                W[Nx_cells + fict + g][j] =
                	W[fict + g][j];
        }
    }

    // ---- Y-границы ----
    for (size_t i = 0; i < Nx_tot; ++i) {
        for (size_t g = 0; g < fict; ++g) {

            // ---- Down ----
            if (y_down_bound == "free")
                W[i][g] = W[i][fict];

            else if (y_down_bound == "wall") {
                W[i][g] = W[i][fict];
                W[i][g][2] *= -1.0;  // меняем v
            }

            else if (y_down_bound == "periodic")
                W[i][g] = W[i][Ny_cells + g];

            // ---- Up ----
            if (y_up_bound == "free")
                W[i][Ny_cells + fict + g] = 
					W[i][Ny_cells + fict - 1];

            else if (y_up_bound == "wall") {
                W[i][Ny_cells + fict + g] =
                    W[i][Ny_cells + fict - 1];
                W[i][Ny_cells + fict + g][2] *= -1.0;
            }

            else if (y_up_bound == "periodic")
                W[i][Ny_cells + fict + g] =
                    W[i][fict + g];
        }
    }
	
	return;

}

