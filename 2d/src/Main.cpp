#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <filesystem>

#include "FileProcessing.h"
#include "Init.h"
#include "BoundCond.h"
#include "RiemannSolver.h"
//#include "Acoustic.h"
#include "GeneralFunctions.h"
//#include "Analis.h"
#include "Limiters.h"
#include "Types.h" // Здесь лежить NEQ - кол-во переменных в состоянии и field


int Nx, Ny;
int step_fo, step_max, bound_case;

double Lx, Ly, t_max, time_fo, x0, gamm, CFL, Q, C1, C2;

std::string x_left_bound, x_right_bound,
			y_up_bound, y_down_bound;

std::string Test, high_order_method, TVD_solver, TVD_limiter;

std::string method, solver, time_method, rec_limiter;
bool Diffusion_flag, Viscous_flag, TVD_flag;

int fict = 3;


void GetDt(const Field& W, 
		   const std::vector<double>& x, 
		   const std::vector<double>& y, 
		   double& dt) {
	
	/*
	double c = 0, dx = Lx, u = 0;
	double c_temp;

	for (int i = fict; i < Nx + fict; i++){
		c_temp = std::sqrt(gamm * std::max(1e-7, W[i][2]) / std::max(1e-7, W[i][0]));
		if (c_temp > c) 
			c = c_temp;
		if (std::abs(W[i][1]) > u)
			u = std::abs(W[i][1]);
		if ((x[i + 1] - x[i]) < dx)
			dx = x[i + 1] - x[i];
	}

	dt = CFL * dx/(c + u);
	return;*/
	double dx = Lx / (Nx - 1);
	double dy = Ly / (Ny - 1);

	double max_lambda_x = 0.0;
    double max_lambda_y = 0.0;

	double rho, u, v, P;
	double c;

	int Nx_cells = Nx - 1;
	int Ny_cells = Ny - 1;

	for (size_t i = fict; i < Nx_cells + fict; i++) {
        for (size_t j = fict; j < Ny_cells + fict; ++j) {

            rho = W[i][j][0];
            u   = W[i][j][1];
            v   = W[i][j][2];
            P   = W[i][j][NEQ - 1];

            c = std::sqrt(gamm * P / rho);

            max_lambda_x = std::max(max_lambda_x,
                                    std::abs(u) + c);

            max_lambda_y = std::max(max_lambda_y,
                                    std::abs(v) + c);
        }
    }

	double dt_x = dx / max_lambda_x;
    double dt_y = dy / max_lambda_y;

	dt = CFL * std::min(dt_x, dt_y);
	return;
}


std::map<std::string, LimiterFunction> get_limiter_map() {
	
    std::map<std::string, LimiterFunction> limiter_map = {
        {"superbee",   &superbeeLim}, 
		{"vanAlbada1", &vanAlbada1}, 
		// Дополнительные лимитеры:
		{"CHARM",      &CHARM},
		{"HCUS",       &HCUS},
		{"HQUICK",     &HQUICK},
		{"Koren",      &Koren},
		{"minmod",     &minmodLim}, // Часто называют просто minmod
		{"MC",         &MC},        // MC (Monotonized Central)
		{"Osher",      &OsherLim},  // Osher and Chakravarthy
		{"ospre",      &ospre},     // Ospre (Oshers-Sweby-P. R. E.)
		{"smart",      &smart},
		{"Sweby",      &Sweby},
		{"UMIST",      &UMIST},
		{"vanAlbada2", &vanAlbada2},
		{"vanLeer",    &vanLeerLim}
    };
    return limiter_map;
}

std::map<std::string, RecLimiterFunction> get_rec_limiter_map() {
	
    std::map<std::string, RecLimiterFunction> limiter_map = {
        {"superbee",   &superbee}, 
		{"minmod",     &minmod}, // Часто называют просто minmod
		{"vanLeer",    &vanLeer}
    };
    return limiter_map;
}

/*
void printProgressBar(double t, double t_max) {
    const int total_cells = 10; // всего 10 ячеек
    int filled_cells = static_cast<int>((t / t_max) * total_cells);

    std::cout << "\r";

    for (int i = 0; i < total_cells; ++i) {
        if (i < filled_cells)
            std::cout << "\033[92m\u25A0\033[0m"; // заполненная ячейка
        else
            std::cout <<"▱"; // пустая ячейка
    }


    int percent = static_cast<int>((t / t_max) * 100);
    std::cout << " " << percent << "%";
    std::cout.flush();
}


*/

int main(int argc, char* argv[]) {

	if (argc < 1) {
        std::cerr << "Используйте: ./test <config_set_name>\n";
        return 1;
    }
    fs::path config_path = argv[1];
    std::string config_name = config_path.stem().string();


	fs::path base_output = fs::path("output") / config_name;

	if (fs::exists(base_output)) {
    	std::cout << "Очищаем папку: " << base_output << std::endl;
    	fs::remove_all(base_output);
	}

	fs::path csv_folder = base_output / "CSV";
	fs::path pics_folder = base_output / "pics";

	fs::create_directories(csv_folder);
	fs::create_directories(pics_folder);

	
	readConfig(config_path.string());


	std::vector<double> xc(Nx + 2 * fict - 1);
	std::vector<double> x(Nx + 2 * fict);	

	std::vector<double> yc(Ny + 2 * fict - 1);
	std::vector<double> y(Ny + 2 * fict);	
	
	Grid(x, xc, y, yc);

	//double rho_L, u_L, P_L;
	//double rho_R, u_R, P_R;
		
	//std::vector<std::vector<double>> W_0(Nx + 2*fict - 1);
	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

	Field W_0(Nx_tot, std::vector<State>(Ny_tot));
	InitValues(W_0, x, y, config_path);
	BoundCond(W_0);

	double t = 0.0, dt = 1.0;

	//SaveFieldToCSV(W_0, x, y, t, (DataFolder / "InitialData.csv").string());
	//std::cout << t_max;

	Field W = W_0;
	Field W_new = W_0;

	int step = 0;
	while (step <= step_max && t <= t_max) {
		GetDt(W, x, y, dt);

		if (t + dt > t_max) 
			dt = t_max - t;
		
		t += dt;
		UpdateArrays(W, 
					W_new,
				    method, 
					high_order_method, 
					solver, 
					TVD_solver,
					&minmodLim,
					&minmod,
					Viscous_flag, 
					time_method, 
					x, y,
					dt);

		BoundCond(W);

		if (step % step_fo == 0)
			SaveFieldToCSV(W, x, y, t,
    			(csv_folder / (std::to_string(step) + "_step.csv")).string());

		
		if (t >= t_max) {
			SaveFieldToCSV(W, x, y, t,
    			(csv_folder / (std::to_string(step) + "_step.csv")).string());
			
			SaveFieldToCSV(W, x, y, t,
    			(csv_folder / "Final.csv").string());
			break;
		}
		
		step++;
	}

/*
	std::vector<double> W_star_ac = {0, 0};
	std::vector<double> W_ac = {0, 0, 0};
	std::vector<std::vector<double>> W_ac1(Nx + 2 * fict -1);
	InitialZeros(W_ac1, 3);
	double e_ac;

	//t_max = t_riemann;	
	

	// I. Решатель Римана
	NewtonForPressure(W_L_riemann, W_R_riemann, W_star_riemann, 1e-6);
	
	std::ofstream file("CSV_files/Riemann.csv");

	file << "x,rho,u,P,e" << std::endl;

	for (int i = fict; i < Nx + fict; i++){
		W_riemann = GetParamsFromChoosingWave(
				W_L_riemann, 
				W_R_riemann, 
				W_star_riemann, 
				xc[i] - x0, 
				t_max);
		e_riemann = W_riemann[2] / (W_riemann[0] * (gamm - 1));
		file << xc[i]  << "," << W_riemann[0] 
			       << "," << W_riemann[1] 
			       << "," << W_riemann[2] 
			       << "," << e_riemann << std::endl;
	}
	file.close();

	// For all methods
	double t = 0.0, rec_time = 0.0;
	double dt_common;
	int step = 0;
*/
	

	/*while (t <= t_max && step <= step_max) {
		
		//GetDt(W_WENO, xc, dt_common);
		GetCommonDt(W_ByMethods, methods, x, dt_common);
		//return 0;
		
		if (t + dt_common > t_max) {
			dt_common = t_max - t;
		}
		t += dt_common;
		
		//printProgressTree(t, t_max);
		printProgressBar(t, t_max);
		LimiterFunction selected_limiter = get_limiter_map().at(TVD_limiter);
		
		for (std::string& method_name : methods) {
			RecLimiterFunction selected_Reclimiter = get_rec_limiter_map().at(RecLimiters_ByMethods[method_name]);
			UpdateArrays(W_ByMethods[method_name], 

				     	 W_new_ByMethods[method_name], 
				     	 W_b_ByMethods[method_name], 
				     	 F_ByMethods[method_name], 
				     	 method_name, 
						 high_order_method,
						 Solver_ByMethods[method_name], 
						 TVD_solver,
						 selected_limiter,
						 selected_Reclimiter,

						 Viscous_flag, 
				     	 TimeIntergation_ByMethods[method_name], x, dt_common);

			BoundCond(W_ByMethods[method_name]);

			if (Diffusion_flag == true) {
				std::vector<std::vector<double>>& current_W = W_ByMethods[method_name];
				std::vector<std::vector<double>> DW(current_W.size());
				InitialZeros(DW, 3);
				std::vector<std::vector<double>> NW(current_W.size());
				InitialZeros(NW, 3);
			
				DOperator(DW, current_W, Q);
				NOperator(NW, current_W, Q); 
			
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < current_W.size(); i++) {

						double original = current_W[i][j];
						double diff_part = DW[i][j];
						double anti_part = NW[i][j];
						double res = original + diff_part + anti_part;
						if (j == 0 || j == 2) {
							if (res <= 1e-7) {
								res = original + diff_part;

								if (res <= 1e-7) {
									res = 1e-7;
								}
							}
						}

						current_W[i][j] = res;
					}
				}
				BoundCond(W_ByMethods[method_name]);
			}

		}
		if (step % step_fo == 0){
			for (std::string& method_name : methods) {
				std::string filename = 
					"CSV_files/ActualRes/StepRec/" 
					+ method_name + "/step_"  
					+ std::to_string(step) 
					+ ".csv";
				file.open(filename);
				WriteToCSV(W_ByMethods[method_name], 
					   xc, 
					   t, 
					   file);	
				file.close();
			 	SaveAnalysisData(AnalysisFiles[method_name], 
                    				t, 
 						W_ByMethods[method_name], 
						xc, 
						W_L_riemann, W_R_riemann, W_star_riemann, x0,
						GetParamsFromChoosingWave,
						RatioL);
			}	
		}		

		if (t > rec_time){
			rec_time += time_fo;
			for (std::string& method_name : methods) {
				std::string filename = 
					"CSV_files/ActualRes/TimeRec/" 
					+ method_name + "/step_"  
					+ std::to_string(step) 
					+ ".csv";
				file.open(filename);
				WriteToCSV(W_ByMethods[method_name], xc, t, file);	
				file.close();
			 	SaveAnalysisData(AnalysisFiles[method_name], 
                    			t, 
 								W_ByMethods[method_name], 
								xc, 
								W_L_riemann, W_R_riemann, W_star_riemann, x0,
								GetParamsFromChoosingWave,
								RatioL);
			}	
		}	

		if (t == t_max){

			for (std::string& method_name : methods) {
				std::string filename = 
					"CSV_files/" 
					+ method_name 
					+ ".csv";
				file.open(filename);
				WriteToCSV(W_ByMethods[method_name], 
					   xc, 
					   t, 
					   file);	
				file.close();
			 	SaveAnalysisData(AnalysisFiles[method_name], 
                    				t, 
 						W_ByMethods[method_name], 
						xc, 
						W_L_riemann, W_R_riemann, W_star_riemann, x0,
						GetParamsFromChoosingWave,
						RatioL);			
			}
			break;			
		}
		step++;
	} 

	for (auto& pair : AnalysisFiles) {
        	pair.second.close();
	}
	*/
	std::cout << std::endl << "Завершено успешно." << std::endl;	
	return 0;

}
