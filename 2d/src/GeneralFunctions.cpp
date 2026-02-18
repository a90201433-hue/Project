#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "RiemannSolver.h"
//#include "ApproxSolvers.h"
#include "BoundCond.h"
#include "TransformValues.h"
#include "Init.h"
#include "GeneralFunctions.h"
//#include "ENO.h"
//#include "WENO.h"
#include "Limiters.h"
#include "Types.h"
#include "FileProcessing.h"

#include "FLIC.h"


extern double gamm, Lx, Ly, C1, C2;
extern int Nx, Ny, fict;

extern std::string method, solver, time_method;
extern std::string high_order_method, TVD_solver;
extern bool Viscous_flag, TVD_flag;

//2D, по индексам, только точный решатель
void SolveBoundProblem(Field W, 
		       		   Field& W_b,
		    		   Field& F,
               		   int dir) {
	
	if (solver == "Exact") {
		
        if (dir == 0){	//циклы переставляются местами, чтобы считалось полосой и только потом переходило на следующую 

            for (int j = fict; j < Ny - 1 + fict; j++) {
                for (int i = fict; i < Nx + fict; i++) {
                    NewtonForPressure(W[i - 1][j], W[i][j], W_b[i][j], 1e-6, dir);
                }
            }
            for (int j = fict; j < Ny - 1 + fict; j++) {
                for (int i = fict; i < Nx + fict; i++) {
                	W_b[i][j] = GetParamsFromChoosingWave(W[i - 1][j], W[i][j], W_b[i][j], 0.0, 1.0, dir);
                }
            }
            return;
        }


        if (dir == 1) {

            for (int i = fict; i < Nx - 1 + fict; i++) {
                for (int j = fict; j < Ny + fict; j++) {
                    NewtonForPressure(W[i][j - 1], W[i][j], W_b[i][j], 1e-6,dir);
                }
            }
            for (int i = fict; i < Nx - 1 + fict; i++) {
                for (int j = fict; j < Ny + fict; j++) {
                	W_b[i][j] = GetParamsFromChoosingWave(W[i][j - 1], W[i][j], W_b[i][j], 0.0, 1.0,dir);
                }
            }
            return;
        }
    }
    /*
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
	}*/

}

//1D , левый правый))
/*void SolveBoundProblem(Field W_L, 
		       		   Field W_R,
		       		   Field& W_b,
		       		   Field& F,
		       		   std::string solver){
	
	if (solver == "Exact") {
        for (int dir = 0; dir <= 1; dir++) {
            for (int i = fict; i < N + fict; i++) {
                for (int j = fict; j < N + fict; j++) {
                    NewtonForPressure(W_L[i], W_R[i], W_b[i], 1e-6, dir);
                }
            }
        }
        for (int dir = 0; dir <= 1; dir++) {
            for (int i = fict; i < N + fict; i++) {
                for (int j = fict; j < N + fict; j++) {
                W_b[i] = GetParamsFromChoosingWave(W_L[i], W_R[i], W_b[i], 0.0, 1.0);
                }
            }
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
}*/


//2D , только Годунов
void FindBoundValues(Field W, 
					 Field& W_b,
					 Field& F,
                     std::vector<double> x,
					 std::vector<double> y,
					 double dt,
					 RecLimiterFunction func,
                     int dir) {

	/*std::vector<std::vector<double>> Slope(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_L(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_R(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_L(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_R(N + 2 * fict, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_tilde(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> W_half(N + 2 * fict - 1, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_half(N + 2 * fict - 1, std::vector<double> (3, 0.0));*/
	

	if (method == "Godunov") {
		SolveBoundProblem(W, W_b, F, dir);
		return;
	}

/*
	else if (method == "Kolgan") {
	
		FindSlopes(W, Slope);	
		ReconstructValues(W, Slope, W_L, W_R, func);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);
		return;		
	}
	
	else if (method == "Kolgan2") {

		ConvertWtoU(W, U);
		FindSlopes(U, Slope);
		ReconstructValues(U, Slope, U_L, U_R, func);
		ConvertUtoW(W_L, U_L);
		ConvertUtoW(W_R, U_R);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);

		return;		
	}

	else if (method == "Rodionov") {
		
		FindSlopes(W, Slope);
		ReconstructValues(W, Slope, W_L, W_R, func);
		TimePredictor(W_tilde, W, W_L, W_R, fict, N + fict - 1, x, dt);
		BoundCond(W_tilde);
		
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < N + 2 * fict - 1; i++){
				W_half[i][j] = 0.5 * (W[i][j] + W_tilde[i][j]);	
			}
		}
		
		ReconstructValues(W_half, Slope, W_L, W_R, func);
		SolveBoundProblem(W_L, W_R, W_b, F, solver);
		return;	
	}
	
	else if (method == "Rodionov2") {
		
		ConvertWtoU(W, U);
		FindSlopes(U, Slope);	
		ReconstructValues(U, Slope, U_L, U_R, func);
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
		ReconstructValues(U_half, Slope, U_L, U_R, func);
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
	}*/
}

//2D

void Streams(Field W,
			 Field& F,
    		 int dir) {

    //F.resize(Nx);

    if (dir == 0) {
        for (int i = fict; i < Nx + fict; i++)
        //F[i].resize(Ny);
            for (int j = fict; j < Ny - 1 + fict; j++) {
			
                F[i][j][0] = W[i][j][0] * W[i][j][1]; // rho*u

                F[i][j][1] = W[i][j][0] * std::pow(W[i][j][1], 2) + W[i][j][NEQ - 1]; // rho*u^2 + P

                F[i][j][2] = W[i][j][0] * W[i][j][1] * W[i][j][2]; // rho*u*v

                double E = 0.5 * W[i][j][0] * (std::pow(W[i][j][1], 2) + std::pow(W[i][j][2], 2)) + W[i][j][NEQ - 1] /(gamm - 1.0);
                F[i][j][NEQ - 1] = W[i][j][1] * (E + W[i][j][NEQ - 1]); // u(rho*E + P)
        }
    }


    if (dir == 1) {

        for (int i = fict; i < Nx - 1 + fict; i++)
        //Flux[i].resize(Ny);
            for (int j = fict; j < Ny + fict; j++) {
                F[i][j][0] = W[i][j][0] * W[i][j][2]; // rho*v

                F[i][j][2] = W[i][j][0] * W[i][j][1] * W[i][j][2]; // rho*u*v

                F[i][j][1] = W[i][j][0] * std::pow(W[i][j][2], 2) + W[i][j][NEQ - 1]; // rho*v^2 + P

                double E = 0.5 * W[i][j][0] * (std::pow(W[i][j][1], 2) + std::pow(W[i][j][2], 2)) + W[i][j][NEQ - 1] /(gamm - 1.0); 
                F[i][j][NEQ - 1] = W[i][j][2] * (E + W[i][j][NEQ - 1]); // u(rho*E + P)
        }
    }
}

//1D
void GetFluxes(// закомментил блок WENO и вязкости
	Field W,
	Field& F,
	RecLimiterFunction func,
	std::vector<double> x,
    std::vector<double> y,
	double dt, 
    int dir) {
		
	if (method == "WENO") {/*
		WENOStreams(W, F);
		if (Viscous_flag) {
			double q = 0.0;
			for (int i = fict + 1; i < N + fict - 1; i++) {
				Visc(W[i - 1], W[i], x[i - 1], x[i], q);
				F[i][1] += q;
			}
		}*/
	} else {

		Field W_b(Nx + 2*fict, std::vector<State>(Ny + 2*fict));
		FindBoundValues(W, W_b, F, x, y, dt, func, dir);

		if (solver == "Exact"){
			Streams(W_b, F, dir);
			/*if (Viscous_flag) {
				double q = 0.0;
				for (int i = fict + 1; i < N + fict - 1; i++) {
					Visc(W[i - 1], W[i], x[i - 1], x[i], q);
					F[i][1] += q;
				}
			}*/
		}
	}
}

void Euler(Field& W_new, 
			Field W, 
			LimiterFunction phi,
			RecLimiterFunction func,
			int init_idx_x, 
			int end_idx_x, 
			int init_idx_y, 
			int end_idx_y, 
			std::vector<double> x, 
			std::vector<double> y, 
			double dt) {
				
	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

	Field F(Nx_tot + 1, std::vector<State>(Ny_tot));
	Field G(Nx_tot,     std::vector<State>(Ny_tot + 1));

	if (method == "TVD") {/*                                      пока без TVD
		Field F_low(N + 2 * fict, {0.0, 0.0, 0.0});
		GetFluxes(W, F_low, "Godunov", TVD_solver, func, x, dt, Viscous_flag);

		Field F_high(N + 2 * fict, {0.0, 0.0, 0.0});
		GetFluxes(W, F_high, high_order_method, TVD_solver, func, x, dt, Viscous_flag);

		std::vector<double> r_array(3, 0.0);

		for (int i = fict; i < N + fict; i++) {
			r_array = FindR(W[i - 1], W[i], W[i+1]);
			for (int j = 0; j < 3; j++) {
				F[i][j] = F_low[i][j] - phi(r_array[j])*(F_low[i][j] - F_high[i][j]);
			}
		}*/

	} /*else if (method == "FLIC") {
		static bool swap = false;   // чередование направлений

		Field W_tilde = W;
		Field W_tmp = W;
		Field W_tmp_2 = W;

		double dt_half = 0.5 * dt;
		if (!swap) {
			FLIC_L(W, W_tilde, x, y, dt_half, 0);
			FLIC_E(W_tilde, W_tmp, x, y, dt_half, 0);

			FLIC_L(W_tmp, W_tilde, x, y, dt, 1);
			FLIC_E(W_tilde, W_tmp_2, x, y, dt, 1);

			FLIC_L(W_tmp_2, W_tilde, x, y, dt_half, 0);
			FLIC_E(W_tilde, W_new, x, y, dt_half, 0);

		} else {
			FLIC_L(W, W_tilde, x, y, dt_half, 1);
			FLIC_E(W_tilde, W_tmp, x, y, dt_half, 1);

			FLIC_L(W_tmp, W_tilde, x, y, dt, 0);
			FLIC_E(W_tilde, W_tmp_2, x, y, dt, 0);

			FLIC_L(W_tmp_2, W_tilde, x, y, dt_half, 1);
			FLIC_E(W_tilde, W_new, x, y, dt_half, 1);
		}

		swap = !swap;

		return;
	
	} */ 
	else {
		GetFluxes(W, F, func, x, y, dt, 0);
        GetFluxes(W, G, func, x, y, dt, 1);
	}

	Field U(Nx_tot, std::vector<State>(Ny_tot));
	ConvertWtoU(W, U, 0);

	Field U_new(Nx_tot, std::vector<State>(Ny_tot));

	for (int i = init_idx_x; i < end_idx_x; i++) {
        double dx = x[i + 1] - x[i];

        for (int j = init_idx_y; j < end_idx_y; j++) {
            double dy = y[j + 1] - y[j];
			for(int k = 0; k < NEQ; k++) {
				//U_new[i][j][k] = U[i][j][k] - dt/dx*(F[i + 1][j][k]- F[i][j][k])- dt/dy*(G[i][j + 1][k]- G[i][j][k]);
				//std::cout << F[i + 1][j][k] << "\n";
			}
			W_new[i][j][0] = std::max(1e-7,
							W[i][j][0]
							- dt/dx * (F[i + 1][j][0] - F[i][j][0])
							- dt/dy * (G[i][j + 1][0] - G[i][j][0]));


            double mom_x_new = W[i][j][0] * W[i][j][1]
							- dt/dx * (F[i + 1][j][1] - F[i][j][1])
							- dt/dy * (G[i][j + 1][1] - G[i][j][1]);  

			W_new[i][j][1] = mom_x_new / W_new[i][j][0];
			if (W_new[i][j][1] < 1e-9) W_new[i][j][1] = 0.0;

			double mom_y_new = W[i][j][0] * W[i][j][2]
				- dt/dx * (F[i + 1][j][2] - F[i][j][2])
				- dt/dy * (G[i][j + 1][2] - G[i][j][2]);

			W_new[i][j][2] = mom_y_new / W_new[i][j][0];
			if (W_new[i][j][2] < 1e-9) W_new[i][j][2] = 0.0;
			
			double E_old =
					W[i][j][3] / (gamm - 1.0)
					+ 0.5 * W[i][j][0] *
					(W[i][j][1]*W[i][j][1] + W[i][j][2]*W[i][j][2]);

			double E_new =
					E_old
					- dt/dx * (F[i + 1][j][3] - F[i][j][3])
					- dt/dy * (G[i][j + 1][3] - G[i][j][3]);

			double kinetic_new =
					0.5 * W_new[i][j][0] *
					(W_new[i][j][1]*W_new[i][j][1] +
					W_new[i][j][2]*W_new[i][j][2]);

			double P_new =
					(gamm - 1.0) * (E_new - kinetic_new);

			W_new[i][j][NEQ - 1] = std::max(1e-6, P_new);

		}


    }
	//std::cout << std::endl;
	//ConvertUtoW(W_new, U_new, 0);
	/*for (int i = init_idx_x; i < end_idx_x; i++) {
		for (int j = init_idx_y; j < end_idx_y; j++) {
            
			for(int k = 0; k < NEQ; k++) {

				std::cout << W_new[i][j][k] << ", "<< U_new[i][j][k] <<  "\n";
		
			}
		}
	}*/
	//SaveFieldToCSV(U_new, x, y, dt, "output/config1/U_new.csv");
	//SaveFieldToCSV(U, x, y, dt, "output/config1/U.csv");
	//SaveFieldToCSV(W_new, x, y, dt, "output/config1/W_new.csv");
}


void UpdateArrays(Field& W, 
				  Field W_new,
				  LimiterFunction phi,
				  RecLimiterFunction func,
				  std::vector<double> x, 
				  std::vector<double> y,
				  double dt) {

	//std::vector<std::vector<double>> W_L(N + 2 * fict);
	//std::vector<std::vector<double>> W_R(N + 2 * fict);
	if (method == "FLIC") FLIC(W_new, W, x, y, dt);
	else if (method == "MacCORMACK") {/* //БЕЗ МАККОРМАКА!!!
		MacCORMACK(W, W_new, x, dt, solver);
		return;*/
	} 

	// Для остальных методов
	else if (time_method == "RK3") {// БЕЗ РК3!!!
		/*RK3(W_new, W, method, solver, func, fict, N + fict - 1, x, dt, Viscous_flag);*/
	} else {
		Euler(W_new, W, phi, func, fict, Nx + fict - 1,  fict, Ny + fict - 1, x, y, dt);
	}

    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {

			/*for (int k = 0; k < NEQ; k++) {
				if (W_new[i][j][k] < 1e-6) W_new[i][j][k] = 0;
			}*/
			W[i][j] = W_new[i][j];
        }
    }
}