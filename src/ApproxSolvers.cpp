#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "lib/TransformValues.h"
#include "lib/RiemannSolver.h"


extern double gamm;
extern int N, fict;

//const size_t Central_num = N + 2 * fict - 1;
//const size_t Bound_num = N + 2 * fict;

//using DataArray = std::vector<std::vector<double>>;

std::vector<double> HLL(std::vector<double> W_L, 
	 		std::vector<double> W_R) {
 
	double rho_L = W_L[0];
	double u_L = W_L[1];
	double P_L = W_L[2];

	double rho_R = W_R[0];
	double u_R = W_R[1];
	double P_R = W_R[2];
	
	double a_L = std::sqrt(gamm*P_L/rho_L);
	double a_R = std::sqrt(gamm*P_R/rho_R);

	double S_L = std::min(u_L - a_L, u_R - a_R);
	
	double S_R = std::max(u_L + a_L, u_R + a_R);
	

	std::vector<double> U_L(3, 0.0), U_R(3, 0.0);
	ConvertWtoU(W_L, U_L);
	ConvertWtoU(W_R, U_R);

	std::vector<double> F_L(3, 0.0), F_R(3, 0.0);
	FindFfromW(F_L, W_L);
	FindFfromW(F_R, W_R);
	

	if (S_L >= 0) {
		return F_L;
	}

	if (S_R <= 0) {
		return F_R;
	}


	std::vector<double> F(3, 0.0);
	

	for (int i = 0; i < 3; i++){
		F[i] = (S_R*F_L[i] - S_L*F_R[i] + S_L*S_R*(U_R[i] - U_L[i]))/(S_R - S_L);
	}

	//std::vector<double> W(3, 0.0);	
	//ConvertUtoW(W, U);
	return F;
}

std::vector<double> HLLC(std::vector<double> W_L, 
	 		std::vector<double> W_R){

	double rho_L = W_L[0], u_L = W_L[1], P_L = W_L[2];
	double rho_R = W_R[0], u_R = W_R[1], P_R = W_R[2];
	
	double a_L = std::sqrt(gamm*P_L/rho_L);
	double a_R = std::sqrt(gamm*P_R/rho_R);
	
	double a_star = 0.5*(u_L - u_R) + (a_L - a_R)/(gamm - 1);
	double u_star = 0.5*(a_L + a_R) + 0.25*(gamm - 1)*(u_L - u_R);

	double S_L = std::min(u_L - a_L, u_star - a_star);
	double S_R = std::max(u_R + a_R, u_star + a_star);

	double S_M = (P_R - P_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R))
					/(rho_L*(S_L - u_L) - rho_R*(S_R - u_R));
	
	double E_L = 0.5*std::pow(u_L, 2) + P_L/(rho_L*(gamm - 1));	
	double E_R = 0.5*std::pow(u_R, 2) + P_R/(rho_R*(gamm - 1));

	std::vector<double> U_L(3, 0.0), U_R(3, 0.0);
	ConvertWtoU(W_L, U_L);
	ConvertWtoU(W_R, U_R);
	
	std::vector<double> F_L(3, 0.0), F_R(3, 0.0);
	FindFfromW(F_L, W_L);
	FindFfromW(F_R, W_R);


	if ((S_R - S_L) < 1e-6) {

		std::vector<double> F(3, 0.0);
		for (int i = 0; i < 3; i++){
			F[i] = 0.5*(F_L[i] + F_R[i]);
		}
		std::cout << "blablabla";
		return F;

	}

	if (S_L >= 0.0) {
		return F_L;
	}

	if (S_L < 0.0 && 0.0 <= S_M){
		double rho_star = rho_L * (S_L - u_L)/(S_L - S_M);
		double P_star = P_L + rho_L*(u_L - S_L)*(u_L - S_M);
		double E_star = ((S_L - u_L)*E_L - P_L*u_L + P_star*S_M)/(S_L - S_M);
		//std::vector<double> W_star = {rho_star, S_M, P_star};
		
		std::vector<double> U_star = {rho_star, rho_star*S_M, rho_star*E_star};
		
		//std::cout << S_M << std::endl;
		std::vector<double> F_star_L(3, 0.0);
		for (int i = 0; i < 3; i++){
			F_star_L[i] = F_L[i] + S_L*(U_star[i] - U_L[i]);
		}

		return F_star_L;
	}
	
	if (S_M <= 0.0 && 0.0 < S_R){
		double rho_star = rho_R * (S_R - u_R)/(S_R - S_M);
		double P_star = P_R + rho_R*(u_R - S_R)*(u_R - S_M);
		double E_star = ((S_R - u_R)*E_R - P_R*u_R + P_star*S_M)/(S_R - S_M);
		std::vector<double> U_star = {rho_star, rho_star*S_M, rho_star*E_star};
		//std::vector<double> W_star = {rho_star, S_M, P_star};
		//std::vector<double> U_star(3, 0.0);
		//ConvertWtoU(W_star, U_star);

		std::vector<double> F_star_R(3, 0.0);
		for (int i = 0; i < 3; i++){
			F_star_R[i] = F_R[i] + S_R*(U_star[i] - U_R[i]);
		}

		return F_star_R;
	}

	return F_R;


	

}

std::vector<double> Rusanov(std::vector<double> W_L, 
                            std::vector<double> W_R) {
  
  	double rho_L = W_L[0];
  	double u_L = W_L[1];
  	double P_L = W_L[2];

  	double rho_R = W_R[0];
  	double u_R = W_R[1];
  	double P_R = W_R[2];
  
  	double a_L = std::sqrt(gamm * P_L / rho_L);
  	double a_R = std::sqrt(gamm * P_R / rho_R);
  
  	double omega = 1.0;

  	double S_max_L = std::abs(u_L) + a_L;
  	double S_max_R = std::abs(u_R) + a_R;

  	double S_max = std::max(S_max_L, S_max_R); 
  	double alpha = omega * S_max; 

  	std::vector<double> U_L(3, 0.0), U_R(3, 0.0);
  	ConvertWtoU(W_L, U_L);
  	ConvertWtoU(W_R, U_R);

  	std::vector<double> F_L(3, 0.0), F_R(3, 0.0);
  	FindFfromW(F_L, W_L);
  	FindFfromW(F_R, W_R);
  
  	// F_{L/R} = 0.5 * (F_L + F_R) - 0.5 * alpha * (U_R - U_L)
  	std::vector<double> F(3, 0.0);

  	for (int i = 0; i < 3; i++){
    		F[i] = 0.5 * (F_L[i] + F_R[i]) - 0.5 * alpha * (U_R[i] - U_L[i]);
  	}

  	return F;
}


// Энтропийная коррекция 

double phi_lambda(double lambda, double delta = 1e-6) {
    if (std::abs(lambda) >= delta) {
        return std::abs(lambda);
    } else {
        return (lambda * lambda + delta * delta) / (2.0 * delta);
    }
}

//F{Roe} = 0.5 * [F_L + F_R - SUM_k(phi(lambda_k) * alpha_k * r_k)]

std::vector<double> Roe(std::vector<double> W_L, 
                        std::vector<double> W_R) {

  double rho_L = W_L[0], u_L = W_L[1], P_L = W_L[2];
  double rho_R = W_R[0], u_R = W_R[1], P_R = W_R[2];
  
    double E_L = P_L / (gamm - 1.0) + 0.5 * rho_L * u_L * u_L;
    double H_L = (E_L + P_L) / rho_L; 

    double E_R = P_R / (gamm - 1.0) + 0.5 * rho_R * u_R * u_R;
    double H_R = (E_R + P_R) / rho_R; 

    // Вычисление Roe-средних параметров 
    double rho_sqrt_L = std::sqrt(rho_L);
    double rho_sqrt_R = std::sqrt(rho_R);
    double rho_sum = rho_sqrt_L + rho_sqrt_R;
    
    double u_tilde   = (rho_sqrt_L * u_L + rho_sqrt_R * u_R) / rho_sum;
    double H_tilde   = (rho_sqrt_L * H_L + rho_sqrt_R * H_R) / rho_sum;
    
    // Скорость звука Roe-средняя 
    double a_tilde = std::sqrt((gamm - 1.0) * (H_tilde - 0.5 * u_tilde * u_tilde));

    std::vector<double> lambda_tilde = {
        u_tilde - a_tilde, 
        u_tilde,           
        u_tilde + a_tilde  
    };


    std::vector<double> U_L(3), U_R(3);
    ConvertWtoU(W_L, U_L);
    ConvertWtoU(W_R, U_R);

    std::vector<double> dU(3);
    for (int i = 0; i < 3; ++i) {
        dU[i] = U_R[i] - U_L[i];
    }
    
    std::vector<double> alpha(3);
    double dP = P_R - P_L;
    double du = u_R - u_L;
    
    double rho_tilde_sq = rho_sqrt_L * rho_sqrt_R; 
    
    alpha[1] = dU[0] - dP / (a_tilde * a_tilde); 
    

    double alpha_term = dP / (rho_tilde_sq * a_tilde);
    
    alpha[0] = 0.5 * (alpha_term - du) / a_tilde * rho_tilde_sq * a_tilde; 
    alpha[2] = 0.5 * (alpha_term + du) / a_tilde * rho_tilde_sq * a_tilde; 
    
    double P_term = dP / (a_tilde * a_tilde);
    alpha[0] = 0.5 * (P_term - rho_tilde_sq * du / a_tilde); // alpha_1
    alpha[2] = 0.5 * (P_term + rho_tilde_sq * du / a_tilde); // alpha_3
    alpha[1] = dU[0] - P_term; // alpha_2
    
    std::vector<std::vector<double>> r(3, std::vector<double>(3));
    
    r[0][0] = 1.0;
    r[0][1] = u_tilde - a_tilde;
    r[0][2] = H_tilde - u_tilde * a_tilde;

    r[1][0] = 1.0;
    r[1][1] = u_tilde;
    r[1][2] = 0.5 * u_tilde * u_tilde;

    r[2][0] = 1.0;
    r[2][1] = u_tilde + a_tilde;
    r[2][2] = H_tilde + u_tilde * a_tilde;
    

    std::vector<double> F_L(3, 0.0), F_R(3, 0.0);
    FindFfromW(F_L, W_L);
    FindFfromW(F_R, W_R);

    std::vector<double> F_roe(3, 0.0);
    
    for (int i = 0; i < 3; i++){
        // Диссипативный член 
        double sum_dissipation = 0.0;
        for (int k = 0; k < 3; k++) {
            sum_dissipation += phi_lambda(lambda_tilde[k]) * alpha[k] * r[k][i];
        }

        F_roe[i] = 0.5 * (F_L[i] + F_R[i]) - 0.5 * sum_dissipation;
    }

  return F_roe;
}

void StarValues(std::vector<double> W_L, 
                std::vector<double> W_R, 
                double& rho_13,
                double& u_13,
                double& c_13,
                double& P_13,
                double& rho_23,
                double& u_23,
                double& c_23,
                double& P_23) {
                                                
    double rho_L = W_L[0], u_L = W_L[1], P_L = W_L[2];
  double rho_R = W_R[0], u_R = W_R[1], P_R = W_R[2];
  
  double c_L = std::sqrt(gamm * P_L / W_L[0]);
  double c_R = std::sqrt(gamm * P_R / W_R[0]);
  
  std::vector<double> W_star;
    NewtonForPressure(W_L, W_R, W_star, 1e-6);
    
    P_13 = W_star[0];
    u_13 = W_star[1];
    c_13 = (u_L + 2.0 * c_L / (gamm - 1.0) - u_13) * 0.5 * (gamm - 1.0);
    rho_13 = rho_L * std::pow(P_13 / P_L, 1.0 / gamm);
    
    P_23 = P_13;
    u_23 = u_13;
    c_23 = (u_23 - (u_R - 2.0 * c_R / (gamm - 1.0))) * 0.5 * (gamm - 1.0);
  rho_23 = rho_R * std::pow(P_23 / P_R, 1.0 / gamm);
  
}

void SonicValues(std::vector<double> W_L, 
                 std::vector<double> W_R,
                 double& rho_16,
                 double& u_16,
                 double& c_16,
                 double& P_16,
                 double& rho_56,
                 double& u_56,
                 double& c_56,
                 double& P_56) {
                                                 
    double rho_L = W_L[0], u_L = W_L[1], P_L = W_L[2];
  double rho_R = W_R[0], u_R = W_R[1], P_R = W_R[2];
  
  double c_L = std::sqrt(gamm * P_L / W_L[0]);
  double c_R = std::sqrt(gamm * P_R / W_R[0]);
  
  u_16 = (gamm - 1.0) * (u_L + 2.0 * c_L / (gamm - 1.0)) / (gamm + 1.0);
  c_16 = u_16;
  double S_L = P_L / std::pow(rho_L, gamm);
  rho_16 = std::pow(std::pow(c_16, 2) / (gamm * S_L), 1 / (gamm - 1.0));
  P_16 = P_L * std::pow(c_16 / c_L, 2.0 * gamm / (gamm - 1.0));
  
  
  u_56 = (gamm - 1.0) * (u_R - 2.0 * c_R / (gamm - 1.0)) / (gamm + 1.0);
  c_56 = -u_56;
  double S_R = P_R / std::pow(rho_R, gamm);
  rho_56 = std::pow(std::pow(c_56, 2) / (gamm * S_R), 1 / (gamm - 1.0));
  P_56 = P_R * std::pow(c_56 / c_R, 2.0 * gamm / (gamm - 1.0));
  
}
  
  
std::vector<double> Osher(std::vector<double> W_L, 
                          std::vector<double> W_R) {


  double rho_L = W_L[0], u_L = W_L[1], P_L = W_L[2];
  double rho_R = W_R[0], u_R = W_R[1], P_R = W_R[2];
  
  double rho_13, u_13, c_13, P_13, rho_23, u_23, c_23, P_23,
         rho_16, u_16, c_16, P_16, rho_56, u_56, c_56, P_56;
  
  double c_L = std::sqrt(gamm * P_L / W_L[0]);
  double c_R = std::sqrt(gamm * P_R / W_R[0]);
    
    StarValues(W_L, 
               W_R, 
               rho_13,
               u_13,
               c_13,
               P_13,
               rho_23,
               u_23,
               c_23,
               P_23);
               
    SonicValues(W_L, 
                W_R,
                rho_16,
                u_16,
                c_16,
                P_16,
                rho_56,
                u_56,
                c_56,
                P_56);
                
    std::vector<double> F_L;
    FindFfromW(F_L, {rho_L, u_L, P_L});            
    std::vector<double> F_R;
    FindFfromW(F_R, {rho_R, u_R, P_R});
    
    std::vector<double> F_13;
    FindFfromW(F_13, {rho_13, u_13, P_13});
    std::vector<double> F_23;
    FindFfromW(F_23, {rho_23, u_23, P_23});
    
    std::vector<double> F_16;
    FindFfromW(F_16, {rho_16, u_16, P_16});
    std::vector<double> F_56;
    FindFfromW(F_56, {rho_56, u_56, P_56});
    
    std::vector<double> F_osher = F_R; 
    
    double lam_L = u_L - c_L;
    double lam_13 = u_13 - c_13;
    
    if (lam_L <= 0 && lam_13 <= 0) {
        F_osher = F_osher;
    } else if (lam_L > 0 && lam_13 > 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_13[i] - F_L[i];
        }
    } else if (lam_L > 0 && lam_13 < 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_16[i] - F_L[i];
        }
    } else if (lam_L < 0 && lam_13 > 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_13[i] - F_16[i];
        }
    }
    
    
    lam_13 = u_13;
    double lam_23 = u_13;
    
    if (lam_13 <= 0 && lam_23 <= 0) {
        F_osher = F_osher;
    } else if (lam_13 > 0 && lam_23 > 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_23[i] - F_13[i];
        }
    } 
    
    lam_23 = u_23 + c_23;
    double lam_R = u_R + c_R;
    
    if (lam_23 <= 0 && lam_R <= 0) {
        F_osher = F_osher;
    } else if (lam_23 > 0 && lam_R > 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_R[i] - F_23[i];
        }
    } else if (lam_23 < 0 && lam_R > 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_R[i] - F_56[i];
        }
    } else if (lam_23 > 0 && lam_R < 0) {
        for (int i = 0; i < 3; i++) {
            F_osher[i] += F_56[i] - F_23[i];
        }
    }
    
    
    
  return F_osher;
}

