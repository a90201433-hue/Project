#include <vector>
#include <string>
#include <iostream>
#include "Types.h"
#include "Domain.h"
#include <mpi.h>

extern int Nx, Ny, fict;
extern std::string x_left_bound, x_right_bound,
				   y_up_bound, y_down_bound;


void BoundCond(Field& W, const Domain& dom) {

    size_t Nx_cells = Nx - 1;
    size_t Ny_cells = Ny - 1;

	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

    bool left  = (dom.rx == 0);
    bool right = (dom.rx == dom.px - 1);
    bool down  = (dom.ry == 0);
    bool up    = (dom.ry == dom.py - 1);

    // ---- X границы ----
    for (size_t j = 0; j < Ny_tot; j++) {
        for (size_t g = 0; g < fict; g++) {

            // LEFT
            if (left) {

                if (x_left_bound == "free")
                    W[g][j] = W[fict][j];

                else if (x_left_bound == "wall") {
                    W[g][j] = W[fict][j];
                    W[g][j][1] *= -1.0;
                }

                else if (x_left_bound == "periodic")
                    W[g][j] = W[Nx_cells + g][j];
            }

            // RIGHT
            if (right) {

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
    }

    // ---- Y границы ----
    for (size_t i = 0; i < Nx_tot; i++) {
        for (size_t g = 0; g < fict; g++) {

            // DOWN
            if (down) {

                if (y_down_bound == "free")
                    W[i][g] = W[i][fict];

                else if (y_down_bound == "wall") {
                    W[i][g] = W[i][fict];
                    W[i][g][2] *= -1.0;
                }

                else if (y_down_bound == "periodic")
                    W[i][g] = W[i][Ny_cells + g];
            }

            // UP
            if (up) {

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
    }
}


void ExchangeGhostCells(Field& W, const Domain& dom) {
    int left  = MPI_PROC_NULL;
    int right = MPI_PROC_NULL;
    int down  = MPI_PROC_NULL;
    int up    = MPI_PROC_NULL;

    if (dom.rx > 0)
        left = dom.ry * dom.px + (dom.rx - 1);

    if (dom.rx < dom.px - 1)
        right = dom.ry * dom.px + (dom.rx + 1);

    if (dom.ry > 0)
        down = (dom.ry - 1) * dom.px + dom.rx;

    if (dom.ry < dom.py - 1)
        up = (dom.ry + 1) * dom.px + dom.rx;

    int Nx_cells = dom.Nx_cells;
    int Ny_cells = dom.Ny_cells;

    // Обмен по x
    int count_x = fict * Ny_cells * NEQ;

    std::vector<double> send_left(count_x);
    std::vector<double> recv_left(count_x);

    std::vector<double> send_right(count_x);
    std::vector<double> recv_right(count_x);

    // Левые
    int k = 0;

    for (int g = 0; g < fict; g++)
    for (int j = fict; j < fict + Ny_cells; j++)
    for (int q = 0; q < NEQ; q++)
        send_left[k++] = W[fict + g][j][q];

    // Правые
    k = 0;

    for (int g = 0; g < fict; g++)
    for (int j = fict; j < fict + Ny_cells; j++)
    for (int q = 0; q < NEQ; q++)
        send_right[k++] = W[fict + Nx_cells - 1 + g][j][q];

    // Обмен
    MPI_Sendrecv(send_left.data(), count_x, MPI_DOUBLE,
                 left, 0,
                 recv_right.data(), count_x, MPI_DOUBLE,
                 right, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(send_right.data(), count_x, MPI_DOUBLE,
                 right, 1,
                 recv_left.data(), count_x, MPI_DOUBLE,
                 left, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Распаковывем
    if (right != MPI_PROC_NULL) {

        k = 0;

        for (int g = 0; g < fict; g++)
        for (int j = fict; j < fict + Ny_cells; j++)
        for (int q = 0; q < NEQ; q++)
            W[Nx_cells + fict + g][j][q] = recv_right[k++];
    }

    if (left != MPI_PROC_NULL) {

        k = 0;

        for (int g = 0; g < fict; g++)
        for (int j = fict; j < fict + Ny_cells; j++)
        for (int q = 0; q < NEQ; q++)
            W[g][j][q] = recv_left[k++];
    }

    // Обмен по y
    int count_y = fict * Nx_cells * NEQ;

    std::vector<double> send_down(count_y);
    std::vector<double> recv_down(count_y);

    std::vector<double> send_up(count_y);
    std::vector<double> recv_up(count_y);

  
    k = 0;

    for (int g = 0; g < fict; g++)
    for (int i = fict; i < fict + Nx_cells; i++)
    for (int q = 0; q < NEQ; q++)
        send_down[k++] = W[i][fict + g][q];


    k = 0;

    for (int g = 0; g < fict; g++)
    for (int i = fict; i < fict + Nx_cells; i++)
    for (int q = 0; q < NEQ; q++)
        send_up[k++] = W[i][fict + Ny_cells - 1 + g][q];


    MPI_Sendrecv(send_down.data(), count_y, MPI_DOUBLE,
                 down, 2,
                 recv_up.data(), count_y, MPI_DOUBLE,
                 up, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(send_up.data(), count_y, MPI_DOUBLE,
                 up, 3,
                 recv_down.data(), count_y, MPI_DOUBLE,
                 down, 3,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    if (up != MPI_PROC_NULL) {

        k = 0;

        for (int g = 0; g < fict; g++)
        for (int i = fict; i < fict + Nx_cells; i++)
        for (int q = 0; q < NEQ; q++)
            W[i][Ny_cells + fict + g][q] = recv_up[k++];
    }

    // ---- unpack down ghost ----
    if (down != MPI_PROC_NULL) {

        k = 0;

        for (int g = 0; g < fict; g++)
        for (int i = fict; i < fict + Nx_cells; i++)
        for (int q = 0; q < NEQ; q++)
            W[i][g][q] = recv_down[k++];
    }
}

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

