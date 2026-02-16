#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_

void PressureInitialGuess(std::vector<double> W_L, std::vector<double> W_R, double& P_prev);

void FindValuesOfFunctions(double P, double P_k, double rho_k, double a_k, double& func_k, double& der_func_k);

void NewtonForPressure(std::vector<double> W_L, std::vector<double> W_R, std::vector<double>& W_star, double eps);

std::vector<double> GetParamsFromChoosingWave(std::vector<double> W_L, std::vector<double> W_R, std::vector<double> W_star, double x, double t);

#endif
