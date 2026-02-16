#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include "ParseTOML.h"
#include "Types.h"

extern int Nx, Ny;
extern int step_fo, step_max, bound_case;

extern double Lx, Ly, t_max, time_fo, x0, gamm, CFL, Q, C1, C2;

extern std::string x_left_bound, x_right_bound,
				   y_up_bound, y_down_bound;

extern std::string Test, high_order_method, TVD_solver, TVD_limiter;
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

	Nx = scheme["N_x"].number;
	Ny = scheme["N_y"].number;
	Lx = scheme["L_x"].number;
	Ly = scheme["L_y"].number;

	CFL = scheme["CFL"].number;
	x_left_bound = scheme["x_left_bound"].str;
	x_right_bound = scheme["x_right_bound"].str;
	y_up_bound = scheme["y_up_bound"].str;
	y_down_bound = scheme["y_down_bound"].str;

	step_fo = toml.root["recording"].table["step_fo"].number;
	time_fo = toml.root["recording"].table["time_fo"].number;
	step_max = toml.root["recording"].table["step_max"].number;

}

void Grid(std::vector<double>& x, std::vector<double>& xc,
		  std::vector<double>& y, std::vector<double>& yc) {
	
	// Шаг сетки
	double dx = Lx / (Nx - 1);
	double dy = Ly / (Ny - 1);

	// Заполняем массивы координат (нужны ли нам центральные?)
	for (int i = 0; i < Nx + 2*fict; i++) 
		x[i] = (i - fict)*dx;
	for (int i = 0; i < Nx + 2*fict - 1; i++) 
		xc[i] = x[i] + 0.5*dx;
	

	for (int i = 0; i < Ny + 2*fict; i++) 
		y[i] = (i - fict)*dy;
	for (int i = 0; i < Ny + 2*fict - 1; i++) 
		yc[i] = y[i] + 0.5*dy;

}

void InitialZeros(std::vector<std::vector<double>> &W_zeros, int intern_size) {
	for (int i = 0; i < (int)W_zeros.size(); i++) {
		W_zeros[i].resize(intern_size, 0.0f);
	}
}


void InitValues(Field& W, 
				const std::vector<double>& x, 
				const std::vector<double>& y,
				const std::string& config_path) {

	SimpleToml config, test;
	
	config.load(config_path);
	Test = config.root["simulation"].table["Test"].str;
	
	// Не тест Сода - тут будет функция заполнения
	if (Test == "custom") {
		std::cerr << "Кастомный тест пока не настроен!" << std::endl;
        return;
	} 

	if (!test.load("tests.toml")) {
        std::cerr << "Нет файла с тестами" << std::endl;
        return;
    }
	
	// Тест Сода вдоль оси x (константа по y)
	auto& values = test.root[Test].table;
	//double dx = Lx / (Nx - 1);
	//double dy = Ly / (Ny - 1);
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

	for (size_t i = fict; i < Nx_tot - fict; i++) {
		for (size_t j = fict; j < Ny_tot - fict; j++) {
			double xc = cell_center(i);
			if (xc < x0)
				W[i][j] = {rho_L, u_L, 0.0, P_L};
			else
				W[i][j] = {rho_R, u_R, 0.0, P_R};
		}
	}
}


