#include<vector>
#include <cmath>
#include <iostream>

extern double gamm;

void AcoustPStarAndUStar(std::vector<double> W_L, 
			 std::vector<double> W_R, 
			 std::vector<double>& W_ac_star) {

	double rho_L = W_L[0];
  	double u_L = W_L[1];
  	double P_L = W_L[2];

  	double rho_R = W_R[0];
  	double u_R = W_R[1];
  	double P_R = W_R[2];
	
	/*
  	double a_L = std::sqrt(gamm * P_L / rho_L);
  	double a_R = std::sqrt(gamm * P_R / rho_R);
	*/
  	if (P_L < 1e-7 || P_R < 1e-7) {
    		W_ac_star[0] = 0.0; //P_star
    		W_ac_star[1] = 0.0; //u_star
    		return;
  	}
/*
  	W_ac_star[0] = (rho_R * rho_L * a_R * a_L * ((u_L - u_R) + P_R / (rho_R * a_R) + P_L / (rho_L * a_L))) / (rho_L * a_L + rho_R * a_R);

  	W_ac_star[1] = (u_L * rho_L * a_L + u_R * rho_R * a_R + (P_L - P_R)) / (rho_L * a_L + rho_R * a_R);
	*/
	
	double P_0 = 0.5*(P_R + P_L);
	double rho_0 = 0.5*(rho_R + rho_L);
	double c = std::sqrt(gamm * P_0 / rho_0);
	double H = rho_0*c;

	W_ac_star[0] = 0.5*(P_L + P_R) - 0.5 * H *(u_R - u_L);

  	W_ac_star[1] = 0.5*(u_L + u_R) - 0.5/H*(P_R - P_L);


	//std::cout << W_ac_star[0] << " " << W_ac_star[1] << std::endl;
	return;
}

std::vector<double> ConfigurationAcousticWaves(std::vector<double> W_L, 
					      std::vector<double> W_R, 
					      std::vector<double> W_ac_star, 
					      
					      double x, double t) {

	std::vector<double> W_acoust;
	double rho_L = W_L[0];
  	double u_L = W_L[1];
  	double P_L = W_L[2];

  	double rho_R = W_R[0];
  	double u_R = W_R[1];
  	double P_R = W_R[2];

 	double a_L = std::sqrt(gamm * P_L / rho_L);
 	double a_R = std::sqrt(gamm * P_R / rho_R);
	double P_0 = 0.5*(P_L + P_R);
	double rho_0 = 0.5*(rho_L + rho_R);
	double a = std::sqrt(gamm*P_0/rho_0);
	//std::cout << a_L << " " << a_R << std::endl;
  	//double c_L = std::sqrt(std::abs((P_R - P_L) / (rho_R - rho_L)));

  	if (x / t <= -a) {
		//std::cout << "Before -ct" << std::endl;
    		W_acoust.push_back(rho_L);
    		W_acoust.push_back(u_L);
    		W_acoust.push_back(P_L);

  	} else if ((-a < x / t) && (x / t <= W_ac_star[1])) {
		//std::cout << "After -ct, before u*" << std::endl;
    		W_acoust.push_back(std::max(0.0d, rho_L+(W_ac_star[0] - P_L)/std::pow(a, 2))); //c* = (gamm*P/rho*)^0.5 => rho* = gamm*P/c*^2
    		W_acoust.push_back(W_ac_star[1]); //W_ac_star[1]=u*
    		W_acoust.push_back(std::max(0.0d, W_ac_star[0])); //W_ac_star[0]=P*
  	} else if ((a > x / t) && (x / t > W_ac_star[1])) {
		//std::cout << "After u*, before ct" << std::endl;
    		W_acoust.push_back(std::max(0.0d, rho_R+(W_ac_star[0] - P_R)/std::pow(a, 2))); //c* = (gamm*P/rho*)^0.5 => rho* = gamm*P/c*^2
    		W_acoust.push_back(W_ac_star[1]); //W_ac_star[1]=u*
    		W_acoust.push_back(std::max(0.0d, W_ac_star[0])); //W_ac_star[0]=P*
  	} else if (x / t >= a) {
		//std::cout << "After ct" << std::endl;
    		W_acoust.push_back(rho_R);
    		W_acoust.push_back(u_R);
    		W_acoust.push_back(P_R);
  	}
	
	return W_acoust;
}
