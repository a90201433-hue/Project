#ifndef _APPROX_SOLVERS_H_
#define _APPROX_SOLVERS_H_

std::vector<double> HLL(std::vector<double> W_L, 
	 		std::vector<double> W_R);

std::vector<double> HLLC(std::vector<double> W_L, 
	 		std::vector<double> W_R);

std::vector<double> Rusanov(std::vector<double> W_L, 
                            std::vector<double> W_R);

double phi_lambda(double lambda, double delta = 1e-6);

std::vector<double> Roe(std::vector<double> W_L, 
                        std::vector<double> W_R);

void StarValues(std::vector<double> W_L, 
                std::vector<double> W_R, 
                double& rho_13,
                double& u_13,
                double& c_13,
                double& P_13,
                double& rho_23,
                double& u_23,
                double& c_23,
                double& P_23);

void SonicValues(std::vector<double> W_L, 
                 std::vector<double> W_R,
                 double& rho_16,
                 double& u_16,
                 double& c_16,
                 double& P_16,
                 double& rho_56,
                 double& u_56,
                 double& c_56,
                 double& P_56);

std::vector<double> Osher(std::vector<double> W_L, 
                          std::vector<double> W_R);


#endif
