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


std::vector<double> Osher(std::vector<double> W_L, 
                          std::vector<double> W_R);


#endif
