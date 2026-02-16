#ifndef _INIT_H_
#define _INIT_H_

void readConfig();

void Grid(double dx, std::vector<double>& x, std::vector<double>& xc);

void InitialZeros(std::vector<std::vector<double>> &W_zeros, int intern_size);

void InitValues(std::string Soda, 
	double& rho_L, double& u_L, double& P_L,
	double& rho_R, double& u_R, double& P_R, 
	double& t_max, double& x0, std::vector<double> xc, 
	std::vector<std::vector<double>>& W);

#endif
