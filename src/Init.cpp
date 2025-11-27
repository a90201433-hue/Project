#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include "lib/ParseTOML.h"


extern int N, fo, fict, step_max, bound_case;
extern double L, t_max, x0, gamm, CFL;
extern std::string x_left_bound, x_right_bound, Soda;
extern std::vector<std::string> methods;

void readConfig() {

	SimpleToml toml;
	
	if (!toml.load("config.toml")) {
        	std::cerr << "Нет config-файла" << std::endl;
        	return;
    	}
	
	
	gamm = toml.root["simulation"].table["gamma"].number;
	Soda = toml.root["simulation"].table["Soda_test"].str;

	auto& scheme = toml.root["scheme"].table;
	
	//auto& method_names = scheme["methods"];

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
			std::cout << "\nНе найдены виды схем. Используется схема Годунова." << std::endl;	
        		methods.push_back("Godunov");
		}
	}		
	

	N = scheme["N_x"].number;
	L = scheme["L_x"].number;
	CFL = scheme["CFL"].number;
	x_left_bound = scheme["x_left_bound"].str;
	x_right_bound = scheme["x_right_bound"].str;

	fo = toml.root["recording"].table["fo"].number;
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


