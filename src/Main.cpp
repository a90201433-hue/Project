#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include "lib/FileProcessing.h"
#include "lib/Init.h"
#include "lib/BoundCond.h"
#include "lib/RiemannSolver.h"
#include "lib/Acoustic.h"
#include "lib/GeneralFunctions.h"
#include "lib/Analis.h"
#include "lib/Limiters.h"

using DataArray = std::vector<std::vector<double>>;

int N, fo, step_max, bound_case;
double L, t_max, x0, gamm, CFL, Q, C1, C2;
std::string x_left_bound, x_right_bound, Soda, high_order_method, TVD_solver, TVD_limiter;
std::vector<std::string> methods, solvers;
bool Diffusion_flag, Viscous_flag, TVD_flag;

int fict = 3;

void GetDt(std::vector<std::vector<double>> W, std::vector<double> x, double& dt){
	double c = 0, dx = L, u = 0;
	double c_temp;

	for (int i = fict; i < N + fict; i++){
		c_temp = std::sqrt(gamm * std::max(1e-7, W[i][2]) / std::max(1e-7, W[i][0]));
		if (c_temp > c) 
			c = c_temp;
		if (std::abs(W[i][1]) > u)
			u = std::abs(W[i][1]);
		if ((x[i + 1] - x[i]) < dx)
			dx = x[i + 1] - x[i];
	}

	dt = CFL * dx/(c + u);
	return;
}


void GetCommonDt(std::map<std::string, DataArray> W, std::vector<std::string> methods, std::vector<double> x, double& dt_common) {
	
	std::vector<double> times;
	double temp_dt = 0.0;
	for (const std::string& method_name : methods) {
		GetDt(W[method_name], x, temp_dt);
		times.push_back(temp_dt);
	}

	dt_common = *std::min_element(times.begin(), times.end());
	return;

}

std::map<std::string, LimiterFunction> get_limiter_map() {
    std::map<std::string, LimiterFunction> limiter_map = {
        {"superbee",   &superbee}, 
		{"vanAlbada1", &vanAlbada1}, 
		// Дополнительные лимитеры, которые вы предоставили:
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
		{"vanLeer",    &vanLeer}
    };
    return limiter_map;
}

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

#include <iostream>
#include <random>
#include <vector>
#include <algorithm>

void printProgressTree(double t, double t_max) {
    const int H = 6;
    const int trunk_h = 2;
    const int trunk_w = 3;
    const int total_rows = H + trunk_h + 1;
    static bool first_call = true;

    // Задаём количество лампочек внутри функции
    const int num_lamps = 10; // <-- здесь можно менять число лампочек

    // Сколько всего листьев
    int total_leaves = H * H;

    // --- вычисляем прогресс ---
    double p = t / t_max;
    if (p < 0) p = 0;
    if (p > 1) p = 1;
    int percent = static_cast<int>(p * 100);
    int filled = static_cast<int>(p * total_leaves);

    // --- фаза гирлянды (0..10) ---
    int phase = percent / 10;

    // --- генератор случайных чисел ---
    static std::mt19937 rng(std::random_device{}());
    auto randColor = [&](int x) {
        switch (x % 5) {
            case 0: return "\033[31m"; // красный
            case 1: return "\033[33m"; // жёлтый
            case 2: return "\033[36m"; // голубой
            case 3: return "\033[35m"; // фиолетовый
            default: return "\033[91m"; // ярко-красный
        }
    };

    // --- храним гирлянду ---
    static std::vector<bool> lamps_mask(total_leaves, false);
    static std::vector<int> lamps_color(total_leaves, 0);
    static int current_phase = -1;

    // если перешли в новую фазу — перегенерируем гирлянду!
    if (phase != current_phase) {
        current_phase = phase;

        // сначала очищаем маску
        std::fill(lamps_mask.begin(), lamps_mask.end(), false);

        // выбираем случайные позиции для лампочек
        std::vector<int> indices(total_leaves);
        for (int i = 0; i < total_leaves; ++i) indices[i] = i;
        std::shuffle(indices.begin(), indices.end(), rng);

        for (int i = 0; i < num_lamps && i < total_leaves; ++i) {
            lamps_mask[indices[i]] = true;
        }

        // генерируем случайные цвета для всех листьев
        for (int i = 0; i < total_leaves; ++i) {
            lamps_color[i] = rng();
        }
    }

    // --- на первом вызове создаём пустое пространство ---
    if (first_call) {
        for (int i = 0; i < total_rows; ++i)
            std::cout << "\n";
        first_call = false;
    }

    // --- поднимаем курсор ---
    std::cout << "\033[" << total_rows << "A";

    int max_width = 2 * (H - 1) + 1;
    int counter = 0;

    // --- Листья ---
    for (int i = 0; i < H; ++i) {
        int width = 2 * i + 1;
        int pad = (max_width - width) / 2;
        std::cout << std::string(pad, ' ');
        for (int j = 0; j < width; ++j) {
            bool isLamp = lamps_mask[counter];
            int color = lamps_color[counter];
            if (counter < filled) {
                if (isLamp) std::cout << randColor(color) << "\033[42m" << "●" << "\033[0m";
                else std::cout << "\033[32m█\033[0m";
            } else {
                if (isLamp) std::cout << randColor(color) << "○" << "\033[0m";
                else std::cout << "░";
            }
            counter++;
        }
        std::cout << "\n";
    }

    // --- Ствол ---
    int pad = (max_width - trunk_w) / 2;
    for (int i = 0; i < trunk_h; ++i) {
        std::cout << std::string(pad, ' ');
        std::cout << "\033[38;5;130m";
        for (int j = 0; j < trunk_w; ++j) std::cout << "█";
        std::cout << "\033[0m\n";
    }

    // --- Процент ---
    std::cout << "Прогресс: " << percent << "%\n";
    std::cout.flush();
}

void InitMaps(std::vector<std::string> methods, std::vector<std::string> solvers, DataArray W_0,
	      	  std::map<std::string, std::string>& Solver_Map,
	      	  std::map<std::string, DataArray>& W_map,
      	      std::map<std::string, DataArray>& W_new_map,
	      	  std::map<std::string, DataArray>& W_b_map,
	      	  std::map<std::string, DataArray>& F_map) {

	const size_t Central_num = N + 2 * fict - 1;
    	const size_t Bound_num = N + 2 * fict;
	
	size_t i = 0;
	for (const std::string& method_name : methods) {
        
        	std::cout << "Инициализация массивов для метода: " << method_name << std::endl;
		
		Solver_Map[method_name] = solvers[i];
		i++;

		W_map[method_name] = W_0;
        
        	// --- 2. W_new ---
        	// Создаем элемент в карте и сразу задаем размер
       	 	W_new_map[method_name].resize(Central_num);
        	InitialZeros(W_new_map[method_name], 3);
        
        	// --- 3. W_b ---
        	W_b_map[method_name].resize(Bound_num);
        	InitialZeros(W_b_map[method_name], 3);
        
        	// --- 4. F ---
        	F_map[method_name].resize(Bound_num);
        	InitialZeros(F_map[method_name], 2);
		
    }	
}

int main() {
	MoveToChache("CSV_files/ActualRes", "CSV_files/ChacheRes");
	ClearDirectoryContents("CSV_files/ActualRes");
	
	readConfig();
	std::cout << TVD_solver << std::endl;
	double dx = L / (N - 1);
	std::vector<double> xc(N + 2 * fict - 1);
	std::vector<double> x(N + 2 * fict);	
	Grid(dx, x, xc);

	double rho_L, u_L, P_L;
	double rho_R, u_R, P_R;
		
	std::vector<std::vector<double>> W_0(N + 2*fict - 1);
	InitValues(Soda, rho_L, u_L, P_L, rho_R, u_R, P_R, t_max, x0, xc, W_0);

	for (auto& method: methods){	
		CreateDir("CSV_files/ActualRes", method);

	}
	
	std::map<std::string, DataArray> W_ByMethods;
	std::map<std::string, DataArray> W_new_ByMethods;
	std::map<std::string, DataArray> W_b_ByMethods;
	std::map<std::string, DataArray> F_ByMethods;
	std::map<std::string, std::string> Solver_ByMethods;	
	InitMaps(methods, solvers, W_0, Solver_ByMethods, 
		 W_ByMethods, W_new_ByMethods, W_b_ByMethods, F_ByMethods);

	
	
	std::vector<double> W_L_riemann = {rho_L, u_L, P_L};
	std::vector<double> W_R_riemann = {rho_R, u_R, P_R};
	std::vector<double> W_star_riemann = {0, 0};
	std::vector<double> W_riemann = {0, 0, 0};
	double e_riemann;

	std::map<std::string, std::ofstream> AnalysisFiles;
    	for (const auto& method : methods) {
        	std::string filename = "CSV_files/Analis_" + method + ".csv";
        	AnalysisFiles[method].open(filename);
        	// Пишем заголовок CSV
        	AnalysisFiles[method] << "time,Ratio_rho,Ratio_v,Ratio_P" << std::endl;
    	}

/*
	std::vector<double> W_star_ac = {0, 0};
	std::vector<double> W_ac = {0, 0, 0};
	std::vector<std::vector<double>> W_ac1(N + 2 * fict -1);
	InitialZeros(W_ac1, 3);
	double e_ac;
*/
	//t_max = t_riemann;	
	

	// I. Решатель Римана
	NewtonForPressure(W_L_riemann, W_R_riemann, W_star_riemann, 1e-6);
	
	std::ofstream file("CSV_files/Riemann.csv");

	file << "x,rho,u,P,e" << std::endl;

	for (int i = fict; i < N + fict; i++){
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
	double t = 0.0;
	double dt_common;
	int step = 0;

	while (t <= t_max) {
		
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
			
			UpdateArrays(W_ByMethods[method_name], 

				     	 W_new_ByMethods[method_name], 
				     	 W_b_ByMethods[method_name], 
				     	 F_ByMethods[method_name], 
				     	 method_name, 
						 high_order_method,
						 Solver_ByMethods[method_name], 
						 TVD_solver,
						 selected_limiter,
						 Viscous_flag,
				     	 "Euler", x, dt_common);
				

				    


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
		if (step % fo == 0){
			for (std::string& method_name : methods) {
				std::string filename = 
					"CSV_files/ActualRes/" 
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

	std::cout << std::endl << "Завершено успешно." << std::endl;	
	return 0;

}
