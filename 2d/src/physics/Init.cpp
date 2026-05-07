#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <cmath>

#include "ParseTOML.h"
#include "Types.h"

extern int Nx, Ny;
extern int Nx_glob, Ny_glob;
extern int step_fo, step_max, bound_case;

extern double Lx, Ly, t_max, time_fo, x0, gamm, gamm1, CFL, Q, C1, C2,
				T_init, R_gas, M, P_min, E_act, Z_freq, VISC, MINWT, GASW, MINGRHO;

extern std::string x_left_bound, x_right_bound,
				   y_up_bound, y_down_bound;

extern std::string high_order_method, TVD_solver, TVD_limiter;
extern std::string method, solver, time_method, rec_limiter;
extern bool Diffusion_flag, Viscous_flag, TVD_flag;
extern int fict;

void readConfig(const std::string& config_path) {

	SimpleToml toml;
	
	if (!toml.load(config_path)) {
        std::cerr << "Нет config-файла" << std::endl;
        return;
    }
	
	gamm = toml.root["simulation"].table["gamma"].number;
	gamm1 = toml.root["simulation"].table["gamma1"].number;
	auto& scheme = toml.root["scheme"].table;
	
	method = scheme["method"].str;

	rec_limiter = scheme["rec_limiter"].str;
	solver = scheme["solver"].str;
	time_method = scheme["time_integration_method"].str;



	Diffusion_flag = (scheme["Diffusion"].str == "On") ? true : false;
	Q = (Diffusion_flag == true) ? scheme["Q"].number : 0.0;

	Viscous_flag = (scheme["Viscous"].str == "On") ? true : false;
	C1 = (Viscous_flag == true) ? scheme["C1"].number : 0.0;
	C2 = (Viscous_flag == true) ? scheme["C2"].number : 0.0;
	
	TVD_flag = (scheme["TVD"].str == "On") ? true : false;
	high_order_method = scheme["High_order_method"].str;
	TVD_solver = scheme["High_order_method"].str;
	TVD_limiter = scheme["High_order_method"].str;

	Nx_glob = scheme["N_x"].number;
	Ny_glob = scheme["N_y"].number;
	Lx = scheme["L_x"].number;
	Ly = scheme["L_y"].number;

	CFL = scheme["CFL"].number;
	x_left_bound = scheme["x_left_bound"].str;
	x_right_bound = scheme["x_right_bound"].str;
	y_up_bound = scheme["y_up_bound"].str;
	y_down_bound = scheme["y_down_bound"].str;

	T_init = scheme["T_init"].number;
	R_gas = scheme["R_gas"].number;
	M = scheme["M"].number;
	P_min = scheme["P_min"].number;
	E_act = scheme["E_act"].number;
	Z_freq = scheme["Z_freq"].number;
	VISC = scheme["VISC"].number;
	MINWT = scheme["MINWT"].number;
	GASW = scheme["GASW"].number;
	MINGRHO = scheme["MINGRHO"].number;

	step_fo = toml.root["recording"].table["step_fo"].number;
	time_fo = toml.root["recording"].table["time_fo"].number;
	step_max = toml.root["recording"].table["step_max"].number;

}


void Grid(std::vector<double>& x, std::vector<double>& y,
          int offset_x, int offset_y) {
	
	// Шаг сетки
	double dx = Lx / (Nx_glob - 1);
	double dy = Ly / (Ny_glob - 1);

	// Заполняем массивы координат (нужны ли нам центральные?)
	for (int i = 0; i < Nx + 2*fict; i++) {
        int global_i = offset_x + i - fict;
		x[i] = global_i * dx;
    }

	for (int j = 0; j < Ny + 2*fict; j++) {
        int global_j = offset_y + j - fict;
        y[j] = global_j * dy;
    }

}
void InitValues(Field& W, 
                const std::vector<double>& x, 
                const std::vector<double>& y,
                const std::string& config_path) {

    SimpleToml config, test;
    
    config.load(config_path);
    std::string Test = config.root["simulation"].table["Test"].str;
    std::string direction = config.root["simulation"].table["direction"].str;

    // ============================================================
    // КАСТОМНЫЙ ТЕСТ: Задача Хааса–Штурггеванта
    // ============================================================
    if (Test == "custom") {
		std::cout << "Запуск задачи Хааса–Штурггеванта (пузырь гелия)" << std::endl;
		
		// Параметры ударной волны (M = 1.22)
		double rho0 = 1.29;      // кг/м³, плотность воздуха при н.у.
		double P0 = 101325.0;    // Па, атмосферное давление
		double u0 = 0.0;
		
		double M_shock = 1.22;
		double gamma_air = 1.4;
		
		double P_ratio = 1.0 + (2.0 * gamma_air / (gamma_air + 1.0)) * (M_shock * M_shock - 1.0);
		double P_shock = P0 * P_ratio;
		
		double rho_ratio = ((gamma_air + 1.0) * M_shock * M_shock) / (2.0 + (gamma_air - 1.0) * M_shock * M_shock);
		double rho_shock = rho0 * rho_ratio;
		
		double u_shock = (M_shock * sqrt(gamma_air * P0 / rho0)) * (1.0 - 1.0 / rho_ratio);
		std::cout<<u_shock<<std::endl;
		// Параметры пузыря гелия
		double R_bubble = 0.025;
		double center_x = 0.175;
		double center_y = 0.0445;
		
		double rho_He = 0.138;
		double P_He = P0;
		double u_He = 0.0;
		
		double gamma_He = 1.66666667;
		
		gamm = gamma_air;
		gamm1 = gamma_He;
		
		double shock_position = 0.12;
		
		size_t Nx_tot = Nx + 2*fict - 1;
		size_t Ny_tot = Ny + 2*fict - 1;
		
		auto cell_center_x = [&](size_t i) {
			return 0.5 * (x[i] + x[i + 1]);
		};
		
		auto cell_center_y = [&](size_t j) {
			return 0.5 * (y[j] + y[j + 1]);
		};
		
		for (size_t i = fict; i < Nx_tot - fict; i++) {
			for (size_t j = fict; j < Ny_tot - fict; j++) {
				double xc = cell_center_x(i);
				double yc = cell_center_y(j);
				
				double dist = sqrt((xc - center_x) * (xc - center_x) + 
								(yc - center_y) * (yc - center_y));
				
				bool inside_bubble = (dist < R_bubble);
				bool behind_shock = (xc < shock_position);
				
				if (inside_bubble) {
					// Внутри пузыря гелия
					W[i][j][0] = rho_He;   // плотность
					W[i][j][1] = u_He;     // скорость U
					W[i][j][2] = 0.0;      // скорость V
					W[i][j][3] = 0.0;      // массовая доля (0 = гелий)
					W[i][j][4] = P_He;     // давление
				} 
				else {
					// Снаружи пузыря (воздух)
					if (behind_shock) {
						W[i][j][0] = rho_shock;
						W[i][j][1] = u_shock;
						W[i][j][2] = 0.0;
						W[i][j][3] = 1.0;  // массовая доля (1 = воздух)
						W[i][j][4] = P_shock;
					} 
					else {
						W[i][j][0] = rho0;
						W[i][j][1] = u0;
						W[i][j][2] = 0.0;
						W[i][j][3] = 1.0;  // массовая доля (1 = воздух)
						W[i][j][4] = P0;
					}
				}
			}
		}
		
		t_max = 0.0006;
		
		std::cout << "Пузырь гелия: R = " << R_bubble << " м, центр = (" 
				<< center_x << ", " << center_y << ")" << std::endl;
		std::cout << "t_max = " << t_max * 1e6 << " мкс" << std::endl;
		
		return;
	}

	// ============================================================
	// СТАНДАРТНЫЕ ТЕСТЫ (Sod, Lax, и т.д.)
	// ============================================================
	/*if (!test.load("tests.toml")) {
        std::cerr << "Нет файла с тестами" << std::endl;
        return;
    }
	
	// Тест Сода вдоль оси x (константа по y)
	auto& values = test.root[Test].table;
	double rho_L, u_L, P_L;
	double rho_R, u_R, P_R;

	rho_L = values["rho_L"].number;
	u_L = values["u_L"].number;
	P_L = values["P_L"].number;

	rho_R = values["rho_R"].number;
	u_R = values["u_R"].number;
	P_R = values["P_R"].number;

	x0 = values["x_gap"].number;
	t_max = values["max_t"].number;

	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

	auto cell_center = [&](size_t k) {
    	return 0.5 * (x[k] + x[k+1]);
	};

	if (direction == "x") {
		for (size_t i = fict; i < Nx_tot - fict; i++) {
			for (size_t j = fict; j < Ny_tot - fict; j++) {
				double xc = cell_center(i);
				if (xc < x0) {
					W[i][j] = {rho_L, u_L, 0.0, P_L};
					mass_fraction_global[i][j][0] = 1.0;  // одно вещество
				}
				else {
					W[i][j] = {rho_R, u_R, 0.0, P_R};
					mass_fraction_global[i][j][0] = 1.0;
				}
			}
		}
	}
	else if (direction == "y") {
		for (size_t i = fict; i < Nx_tot - fict; i++) {
			for (size_t j = fict; j < Ny_tot - fict; j++) {
				double yc = cell_center(j);
				if (yc < x0) {
					W[j][i] = {rho_L, 0.0, u_L, P_L};
					mass_fraction_global[j][i][0] = 1.0;
				}
				else {
					W[j][i] = {rho_R, 0.0, u_R, P_R};
					mass_fraction_global[j][i][0] = 1.0;
				}
			}
		}
	}*/
}


