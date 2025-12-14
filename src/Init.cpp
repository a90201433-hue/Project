#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include "lib/ParseTOML.h"


extern int N, step_fo, step_max, bound_case;
extern double L, t_max, time_fo, x0, gamm, CFL, Q, C1, C2;

extern std::string x_left_bound, x_right_bound, Soda, high_order_method, TVD_solver, TVD_limiter;
extern std::vector<std::string> methods, solvers, time_methods;
extern bool Diffusion_flag, Viscous_flag, TVD_flag;
extern int fict;

void readConfig() {

	SimpleToml toml;
	
	if (!toml.load("config.toml")) {
        	std::cerr << "Нет config-файла" << std::endl;
        	return;
    	}
	
	gamm = toml.root["simulation"].table["gamma"].number;
	Soda = toml.root["simulation"].table["Soda_test"].str;

	auto& scheme = toml.root["scheme"].table;
	
	TVD_flag = (scheme["TVD"].str == "On") ? true : false;
	if (TVD_flag) {
		high_order_method = scheme["High_order_method"].str;
		TVD_solver = scheme["TVD_solver"].str;
		methods.push_back("TVD");
		solvers.push_back(TVD_solver);
		time_methods.push_back("Euler");
	}

	TVD_limiter = (TVD_flag == true) ? scheme["TVD_limiter"].str : "superbee";

	if (scheme.count("methods") == 0) {
		std::cout << "\nНе найдены виды схем. Используется схема Годунова." << std::endl;
    		methods = {"Godunov"};
	} else {
    		TomlValue &arr = scheme["methods"];

    		// извлекаем список
    		for (auto &x : arr.list)
        	methods.push_back(x.str);

    		// если список пустой
    		if (methods.empty()){
			std::cout << "\n" << std::endl;	
        		methods.push_back("Godunov");
		}
	}		
	
	if (scheme.count("solvers") == 0) {
		std::cout << "\nНе найдены решатели. Используется приближенный решатель Римана HLL." << std::endl;
    		solvers = {"HLL"};
	} else {
    		TomlValue &arr = scheme["solvers"];

    		// извлекаем список
    		for (auto &x : arr.list)
        	solvers.push_back(x.str);

    		// если список пустой
    		if (methods.empty()){
			std::cout << "\nНе найдены решатели. Используется приближенный решатель Римана HLL." << std::endl;	
        		solvers.push_back("HLL");
		}
	}

	if (methods.size() > solvers.size()){
		size_t n = methods.size(), k = solvers.size();
		for(size_t i = 0; i < (n - k); i++){
			solvers.push_back("HLL");
		}
	}

	if (scheme.count("time_integration_methods") == 0) {
		std::cout << "\nНе найдены методы интегирования. Используется интегрирование по Эйлеру" << std::endl;
    		time_methods = {"Euler"};
	} else {
    		TomlValue &arr = scheme["time_integration_methods"];

    		// извлекаем список
    		for (auto &x : arr.list)
        	time_methods.push_back(x.str);

    		// если список пустой
    		if (time_methods.empty()){
			std::cout << "\nНе найдены методы интегирования. Используется интегрирование по Эйлеру" << std::endl;	
        		time_methods.push_back("Euler");
		}
	}

	if (methods.size() > time_methods.size()){
		size_t n = methods.size(), k = time_methods.size();
		for(size_t i = 0; i < (n - k); i++){
			time_methods.push_back("Euler");
		}
	}


	Diffusion_flag = (scheme["Diffusion"].str == "On") ? true : false;
	Q = (Diffusion_flag == true) ? scheme["Q"].number : 0.0;

	Viscous_flag = (scheme["Viscous"].str == "On") ? true : false;
	C1 = (Viscous_flag == true) ? scheme["C1"].number : 0.0;
	C2 = (Viscous_flag == true) ? scheme["C2"].number : 0.0;
	

	N = scheme["N_x"].number;
	L = scheme["L_x"].number;
	CFL = scheme["CFL"].number;
	x_left_bound = scheme["x_left_bound"].str;
	x_right_bound = scheme["x_right_bound"].str;


	step_fo = toml.root["recording"].table["step_fo"].number;
	time_fo = toml.root["recording"].table["time_fo"].number;
	step_max = toml.root["recording"].table["step_max"].number;

}

void Grid(double dx, std::vector<double>& x, std::vector<double>& xc){
	
	for (int i = 0; i < N + 2*fict; i++){
		x[i] = (i - fict)*dx;
	}

	for (int i = 0; i < N + 2*fict - 1; i++){
		xc[i] = x[i] + 0.5*dx;
	}
}

void InitialZeros(std::vector<std::vector<double>> &W_zeros, int intern_size) {
	for (int i = 0; i < (int)W_zeros.size(); i++) {
		W_zeros[i].resize(intern_size, 0.0f);
	}
}


void InitValues(std::string Soda, 
	double& rho_L, double& u_L, double& P_L,
	double& rho_R, double& u_R, double& P_R, 
	double& t_max, double& x0, std::vector<double> xc, 
	std::vector<std::vector<double>>& W){

	SimpleToml toml;
	
	if (!toml.load("Tests.toml")) {
        	std::cerr << "Нет файла с тестами" << std::endl;
        	return;
    	}
		
	auto& values = toml.root[Soda].table;
		
	rho_L = values["rho_L"].number;
	u_L = values["u_L"].number;
	P_L = values["P_L"].number;

	rho_R = values["rho_R"].number;
	u_R = values["u_R"].number;
	P_R = values["P_R"].number;

	x0 = values["x_gap"].number;
	t_max = values["max_t"].number;
		
	
	for (int i = 0; i < N + 2*fict - 1; i++){
		if (xc[i] < x0)
			W[i] = {rho_L, u_L, P_L};
		
		else
			W[i] = {rho_R, u_R, P_R};
	}
}


