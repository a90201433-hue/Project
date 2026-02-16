#include <iostream>
#include <vector>
#include <cmath>

extern double gamm;

void PressureInitialGuess(std::vector<double> W_L, std::vector<double> W_R, double& P_prev){
	
	double rho_L = W_L[0];
	double u_L = W_L[1];
	double P_L = W_L[2];

	double rho_R = W_R[0];
	double u_R = W_R[1];
	double P_R = W_R[2];

	double Q_user = 2.;
	double P_min = std::min(P_L, P_R);
	double P_max = std::max(P_L, P_R);
	double a_L = std::sqrt(gamm * P_L / rho_L);
	double a_R = std::sqrt(gamm * P_R / rho_R);

	double P_pvrs = std::max(1e-6, 0.5 * (P_L + P_R) - 0.125 * (u_R - u_L) * (rho_L + rho_L) * (a_L + a_R));
	
	double Q = P_max / P_min;

	if ((Q < Q_user) && (P_min < P_pvrs) && (P_pvrs < P_max)) {
		if (P_pvrs < P_min) {
			P_prev = std::pow((a_L + a_R - 0.5*(gamm - 1)*(u_R - u_L))/(a_L/std::pow(P_L, (gamm - 1)/(2*gamm)) + a_R / std::pow(P_R, (gamm - 1) / (2 * gamm))),(2*gamm)/(gamm - 1));
		}
		else {
			double A_L = 2 / ((gamm + 1) * rho_L);
			double B_L = (gamm - 1) * P_L / (gamm + 1);
			double A_R = 2 / ((gamm + 1) * rho_R);
			double B_R = (gamm - 1) * P_R / (gamm + 1);

			double g_L = std::pow(A_L / (P_pvrs + B_L), 0.5);
			double g_R = std::pow(A_R / (P_pvrs + B_R), 0.5);
				
			P_prev = (g_L*P_L + g_R*P_R - (u_R - u_L))/(g_L + g_R);
			P_prev = std::max(1e-6d, P_prev);
		}
	}
	else {
		P_prev = P_pvrs;
	}
}

void FindValuesOfFunctions(double P, double P_k, double rho_k, double a_k, double& func_k, double& der_func_k) {
	
	if (P <= P_k) {
		//std::cout << "P/P_k " <<  P / P_k<< std::endl;
		func_k = 2 * a_k / (gamm - 1) * (std::pow((P / P_k), (gamm - 1) / (2 * gamm)) - 1);
		der_func_k = 1 / (rho_k * a_k) * std::pow((P / P_k), -(gamm + 1) / (2 * gamm));

	}

	else if (P > P_k) {
		double A_k = 2 / ((gamm + 1) * rho_k);
		double B_k = (gamm - 1) * P_k / (gamm + 1);
		func_k = (P - P_k) * std::sqrt(A_k / (P + B_k));
		der_func_k = std::sqrt(A_k / (P + B_k)) * (1 - (P - P_k) / (2 * (P + B_k)));
	}
	
}

void NewtonForPressure(std::vector<double> W_L, std::vector<double> W_R, std::vector<double>& W_star, double eps) {	
/*	
	if (W_L == W_R){
		W_star = {W_L[2], W_R[1]};
		return;
	}
*/
	double rho_L = W_L[0];
	double u_L = W_L[1];
	double P_L = W_L[2];

	double rho_R = W_R[0];
	double u_R = W_R[1];
	double P_R = W_R[2];

	
	if (P_L <= 1e-7 || P_R <= 1e-7){
	  	W_star[0] = 0; //P_star
		W_star[1] = 0; //u_star
		return;
	}


	double func_L, der_func_L, func_R, der_func_R, P_prev, P_new, Phi;

	double a_L = std::sqrt(gamm * P_L / rho_L);
	double a_R = std::sqrt(gamm * P_R / rho_R);

	PressureInitialGuess(W_L, W_R, P_prev);
	
	FindValuesOfFunctions(P_prev, P_L, rho_L, a_L, func_L, der_func_L);
	FindValuesOfFunctions(P_prev, P_R, rho_R, a_R, func_R, der_func_R);
	
	Phi = func_L + func_R + u_R - u_L;
	P_new = std::max(1e-6, P_prev - Phi / (der_func_L + der_func_R));
/*	
	std::cout << "left: " << W_L[0] << " " <<  W_L[1] << " " << W_L[2] << std::endl;  
	std::cout << "right: " << W_R[0] << " " <<  W_R[1] << " " << W_R[2] << std::endl;
	std::cout << P_new << std::endl;*/
	size_t counter = 0;	
	while (std::abs(P_new - P_prev)/(0.5*(std::abs(P_new) + std::abs(P_prev))) > eps) {
		P_prev = P_new;

		FindValuesOfFunctions(P_prev, P_L, rho_L, a_L, func_L, der_func_L);
		
		FindValuesOfFunctions(P_prev, P_R, rho_R, a_R, func_R, der_func_R);

		Phi = func_L + func_R + u_R - u_L;
		P_new = P_prev - Phi / (der_func_L + der_func_R);
		
		if (counter >  20){
			P_new = 0.5 * (P_L + P_R) - 0.125 * (u_R - u_L) * (rho_L + rho_L) * (a_L + a_R);
			std::cout << "Counter works" << std::endl;
			break;
		}
		counter++;
	}
	//std::cout << P_new << std::endl;
	W_star[0] = P_new;
	FindValuesOfFunctions(P_new, P_L, rho_L, a_L, func_L, der_func_L);
	FindValuesOfFunctions(P_new, P_R, rho_R, a_R, func_R, der_func_R);
	W_star[1] = 0.5 * (u_L + u_R + func_R - func_L);
}
// Function returns vector W = (rho, u, P) for definite x and t and corresponding left, right and star values
std::vector<double> GetParamsFromChoosingWave(std::vector<double> W_L, std::vector<double> W_R, std::vector<double> W_star, double x, double t){
    
	double rho_L = W_L[0];
	double u_L = W_L[1];
	double P_L = W_L[2];

	double rho_R = W_R[0];
	double u_R = W_R[1];
	double P_R = W_R[2];
	
	double P_star = W_star[0];
	double u_star = W_star[1];
	

    	double a_L = std::sqrt(gamm * P_L / rho_L);
    	double a_star_L = a_L * std::pow((P_star / P_L), (gamm - 1) / (2 * gamm));
    	double S_star_L = u_L + 2 * a_L / (gamm - 1);

    	double u0 = u_L + 2 * a_L / (gamm - 1);

    	double a_R = std::sqrt(gamm * P_R / rho_R);
    	double a_star_R = a_R * std::pow((P_star / P_R), (gamm - 1) / (2 * gamm));
    	double S_star_R = u_R - 2 * a_R / (gamm - 1);

    	std::vector<double> W(0);
/*
	std::cout << std::endl;
	std::cout << "Next bound" << std::endl;	
	std::cout << "rho_L = " << rho_L << ", " << "u_L = " << u_L << ", " << "P_L = " << P_L << std::endl;
	std::cout << "rho_R = " << rho_R << ", " << "u_R = " << u_R << ", " << "P_R = " << P_R << std::endl;
	std::cout << "P_star = " << P_star << ", " << "u_star = " << u_star << std::endl;
*/
    	// With vacuum
    	if (P_L > 1e-7 && P_R <= 1.1e-7) // Right vacuum
    	{
		//std::cout << "Right vacuum" << std::endl;
		if ((x / t) <= u_L - a_L) // To the left of head rarefaction
        	{
			
			//std::cout << "Vacuum, right area, after rare 1" << std::endl;
            		W.push_back(rho_L);
            		W.push_back(u_L);
            		W.push_back(P_L);
            		//std::cout << "Vacuum, right area, after rare 2" << std::endl;
        	}
        	else if ((u_L - a_L < (x / t)) && ((x / t) < S_star_L)) // Between head rarefaction and vacuum front (fan)
        	{
			//std::cout << "Vacuum, right area, fan 1" << std::endl;
            		W.push_back(rho_L * std::pow(2 / (gamm + 1) + ((gamm - 1) / ((gamm + 1) * a_L)) * (u_L - (x / t)), 2 / (gamm - 1)));
            		W.push_back((2 / (gamm + 1)) * (a_L + ((gamm - 1) / 2) * u_L + x / t));
            		W.push_back(P_L * std::pow(2 / (gamm + 1) + ((gamm - 1) / ((gamm + 1) * a_L)) * (u_L - (x / t)), (2 * gamm) / (gamm - 1)));
           		//std::cout << "Vacuum, right area, fan 2" << std::endl;
        	}
        	else if ((x / t) >= S_star_L) // To the right of vacuum front
        	{
			//std::cout << "Vacuum, right area, after vacuum front 1" << std::endl;
            		W.push_back(0);
            		W.push_back(u0);
            		W.push_back(0);
            		//std::cout << "Vacuum, right area, after vacuum front 2" << std::endl;
        	}
	}

    	else if (P_R > 1e-7 && P_L <= 1.1e-7) // Left vacuum
    	{
		//std::cout << "Left vacuum" << std::endl;
        	if ((x / t) <= S_star_R) // To the left of vacuum front
        	{
			//std::cout << "Vacuum, left area, after vacuum front 1" << std::endl;
            		W.push_back(0);
            		W.push_back(u0);
            		W.push_back(0);
            		//std::cout << "Vacuum, left area, after vacuum front 2" << std::endl;
        	}
       		else if ((S_star_R < (x / t)) && ((x / t) < u_R + a_R)) // Between head rarefaction and vacuum front (fan)
        	{
			//std::cout << "Vacuum, left area, fan 1" << std::endl;
            		W.push_back(rho_R * std::pow(2 / (gamm + 1) - ((gamm - 1) / ((gamm + 1) * a_R)) * (u_R - (x / t)), 2 / (gamm - 1)));
            		W.push_back((2 / (gamm + 1)) * (-a_R + ((gamm - 1) / 2) * u_R + x / t));
            		W.push_back(P_R * std::pow(2 / (gamm + 1) - ((gamm - 1) / ((gamm + 1) * a_R)) * (u_R - (x / t)), (2 * gamm) / (gamm - 1)));
            		//std::cout << "Vacuum, left area, fan 2" << std::endl;
        	}
        	else if ((x / t) >= u_R + a_R) // To the right of head rarefaction
        	{
			//std::cout << "Vacuum, left area, after rare 1" << std::endl;
            		W.push_back(rho_R);
            		W.push_back(u_R);
            		W.push_back(P_R);
            		//std::cout << "Vacuum, left area, after rare 2" << std::endl;
        	}
    	}
	else if ((P_L <= 1.1e-7 && P_R <= 1.1e-7))
    	{
		//else if (2 * (a_L + a_R)/(gamm - 1) <= u_R - u_L) // Generation vacuum
		//std::cout << "Generation vacuum" << std::endl;
        	if ((x / t) <= S_star_L)
        	{
            		if ((x / t) <= u_L - a_L) // To the left of head rarefaction
            		{
				//std::cout << "Generation vacuum, left area, after rare 1" << std::endl;
                		W.push_back(rho_L);
                		W.push_back(u_L);
                		W.push_back(P_L);
                		//std::cout << "Generation vacuum, left area, after rare 2" << std::endl;
            		}
            		else if ((u_L - a_L < (x / t)) && ((x / t) < S_star_L)) // Between head rarefaction and vacuum front (fan)
            		{
				//std::cout << "Generation vacuum, left area, fan 1" << std::endl;
               			W.push_back(rho_L * std::pow(2 / (gamm + 1) + ((gamm - 1) / ((gamm + 1) * a_L)) * (u_L - (x / t)), 2 / (gamm - 1)));
                		W.push_back((2 / (gamm + 1)) * (a_L + ((gamm - 1) / 2) * u_L + x / t));
                		W.push_back(P_L * std::pow(2 / (gamm + 1) + ((gamm - 1) / ((gamm + 1) * a_L)) * (u_L - (x / t)), (2 * gamm) / (gamm - 1)));
                		//std::cout << "Generation vacuum, left area, fan 2" << std::endl;
            		}

            		else if ((x / t) >= S_star_L) // To the right of vacuum front
            		{
				//std::cout << "Generation vacuum, left area, after vacuum 1" << std::endl;
                		W.push_back(0);
                		W.push_back(u0);
                		W.push_back(0);
                		//std::cout << "Generation vacuum, left area, after vacuum 2" << std::endl;
           		}
		}
	
		else if ((S_star_L < (x / t)) && ((x / t) < S_star_R)) // Between two vacuum fronts
        	{
			//std::cout << "Generation vacuum, central area, between two vacuum fronts 1" << std::endl;
            		W.push_back(0);
            		W.push_back(u0);
            		W.push_back(0);
            		//std::cout << "Generation vacuum, central area, between two vacuum fronts 2" << std::endl;
        	}

        	else if ((x / t) >= S_star_R)
        	{
            		if ((x / t) <= S_star_R) // To the left of vacuum front
            		{
				//std::cout << "Generation vacuum, right area, after vacuum 1" << std::endl;
                		W.push_back(0);
                		W.push_back(u0);
                		W.push_back(0);
                		//std::cout << "Generation vacuum, right area, after vacuum 2" << std::endl;
            		}
	    		else if ((S_star_R < (x / t)) && ((x / t) < u_R + a_R)) // Between head rarefaction and vacuum front (fan)
            		{
				//std::cout << "Generation vacuum, right area, fan 1" << std::endl;
                		W.push_back(rho_R * std::pow(2 / (gamm + 1) - ((gamm - 1) / ((gamm + 1) * a_R)) * (u_R - (x / t)), 2 / (gamm - 1)));
                		W.push_back((2 / (gamm + 1)) * (-a_R + ((gamm - 1) / 2) * u_R + x / t));
                		W.push_back(P_R * std::pow(2 / (gamm + 1) - ((gamm - 1) / ((gamm + 1) * a_R)) * (u_R - (x / t)), (2 * gamm) / (gamm - 1)));
                		//std::cout << "Generation vacuum, right area, fan 2" << std::endl;;
            		}
            		else if ((x / t) >= u_R + a_R) // To the right of head rarefaction
            		{
				//std::cout << "Generation vacuum, right area, after rare 1" << std::endl;
                		W.push_back(rho_R);
                		W.push_back(u_R);
                		W.push_back(P_R);
                		//std::cout << "Generation vacuum, right area, after rare 2" << std::endl;
            		}
        	}
    	}
	
	// Without vacuum
    	if ((x / t) <= u_star) // Left area 
    	{
        	if (P_star > P_L) // Left shock
        	{
            		double S_L = u_L - a_L * std::sqrt(((gamm + 1) / (2 * gamm)) * (P_star / P_L) + (gamm - 1) / (2 * gamm)); // ksi shock
            		if ((S_L <= (x / t)) && ((x / t) <= u_star)) // Between left shock and contact discontinuity
            		{	//std::cout << std::endl;
				//std::cout << "Not vacuum, left area, between shock and contact discontinuity 1" << std::endl;
                		double rho_star_L_sho = rho_L * ((P_star / P_L + (gamm - 1) / (gamm + 1)) / (((gamm - 1) / (gamm + 1)) * (P_star / P_L) + 1));
                		W.push_back(rho_star_L_sho);
                		W.push_back(u_star);
                		W.push_back(P_star);
                		//std::cout << "Not vacuum, left area, between shock and contact discontinuity 2" << std::endl;
            		}
            		else if ((x / t) < S_L) // To the left of the left shock
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, left area, after shock 1" << std::endl;
                		W.push_back(rho_L);
                		W.push_back(u_L);
                		W.push_back(P_L);
                		//std::cout << "Not vacuum, left area, after shock 2" << std::endl;
            		}
       		}
        	else // Left rarefaction
        	{
            		double S_HL = u_L - a_L; // ksi head
            		double S_TL = u_star - a_star_L; // ksi tail
            		if ((x / t) <= S_HL) // To the left of the head of left rarefaction
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, left area, after rare 1" << std::endl;
                		W.push_back(rho_L);
                		W.push_back(u_L);
                		W.push_back(P_L);
                		//std::cout << "Not vacuum, left area, after rare 2" << std::endl;
            		}
            		else if ((S_HL < (x / t) && (x / t) < S_TL)) // Between head and tail of left rarfaction (fan)
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, left area, fan 1" << std::endl;
                		W.push_back(rho_L * std::pow(2 / (gamm + 1) + ((gamm - 1) / ((gamm + 1) * a_L)) * (u_L - x / t), 2 / (gamm - 1)));
                		W.push_back((2 / (gamm + 1)) * (a_L + ((gamm - 1) / 2) * u_L + x / t));
                		W.push_back(P_L * std::pow(2 / (gamm + 1) + ((gamm - 1) / ((gamm + 1) * a_L)) * (u_L - x / t), (2 * gamm) / (gamm - 1)));
               			//std::cout << "Not vacuum, left area, fan 2" << std::endl;
            		}
            		else if ((S_TL <= (x / t)) && ((x / t) <= u_star)) // Between tail of left rarefsction and contact discontinuity
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, left area, between rare and contact discontinuity 1" << std::endl;
                		double rho_star_L_fan = rho_L * std::pow(P_star / P_L, 1 / gamm);
                		W.push_back(rho_star_L_fan);
                		W.push_back(u_star);
                		W.push_back(P_star);
                		//std::cout << "Not vacuum, left area, between rare and contaact discontinuity 2" << std::endl;
            		}
        	}
    	}

    	else  // Right area 
    	{
		
        	if (P_star > P_R) // Right shock
        	{
			
            		double S_R = u_R + a_R * std::sqrt(((gamm + 1) / (2 * gamm)) * (P_star / P_R) + (gamm - 1) / (2 * gamm)); // ksi shock
            		if ((S_R >= (x / t)) && ((x / t) >= u_star)) // Between right shock and contact discontinuity
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, right area, between rare and contact discontinuity 1" << std::endl;
                		double rho_star_R_sho = rho_R * ((P_star / P_R + (gamm - 1) / (gamm + 1)) / (((gamm - 1) / (gamm + 1)) * (P_star / P_R) + 1));
                		W.push_back(rho_star_R_sho);
                		W.push_back(u_star);
                		W.push_back(P_star);
                		//std::cout << "Not vacuum, right area, between rare and contact discontinuity 2" << std::endl;
            		}
            		else if ((x / t) > S_R) // To the right of the right shock
            		{
				//std::cout << std::endl;
			 	//std::cout << "Not vacuum, right area, after shock 1" << std::endl;
                		W.push_back(rho_R);
                		W.push_back(u_R);
                		W.push_back(P_R);
                		//std::cout << "Not vacuum, right area, after shock 2" << std::endl;
            		}
        	}
		else // Right rarefaction
        	{
            		double S_HR = u_R + a_R; // ksi head
            		double S_TR = u_star + a_star_R; // ksi tail
            		if ((x / t) >= S_HR) // To the right of the head of right rarefaction
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, right area, after rare 1" << std::endl;
                		W.push_back(rho_R);
                		W.push_back(u_R);
                		W.push_back(P_R);
                		//std::cout << "Not vacuum, right area, after rare 2" << std::endl;
            		}
            		else if ((S_HR > (x / t) && (x / t) >= S_TR)) // Between head and tail of right rarfaction (fan)
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, right area, fan 1" << std::endl;
                		W.push_back(rho_R * std::pow(2 / (gamm + 1) - ((gamm - 1) / ((gamm + 1) * a_R)) * (u_R - x / t), 2 / (gamm - 1)));
                		W.push_back((2 / (gamm + 1)) * (-a_R + ((gamm - 1) / 2) * u_R + x / t));
                		W.push_back(P_R * std::pow(2 / (gamm + 1) - ((gamm - 1) / ((gamm + 1) * a_R)) * (u_R - x / t), (2 * gamm) / (gamm - 1)));
                		//std::cout << "Not vacuum, right area, fan 2" << std::endl;
            		}
            		else if ((S_TR > (x / t)) && ((x / t) >= u_star)) // Between tail of right rarefaction and contact discontinuity
            		{
				//std::cout << std::endl;
				//std::cout << "Not vacuum, right area, between rare and contact discontinuity 1" << std::endl;
                		double rho_star_R_fan = rho_R * std::pow(P_star / P_R, 1 / gamm);
                		W.push_back(rho_star_R_fan);
                		W.push_back(u_star);
                		W.push_back(P_star);
                		//std::cout << "Not vacuum, right area, between rare and contact discontinuity 2" << std::endl;
			}
		}
	}
    
	//std::cout << std::endl;
    	//std::cout << "After GetParams: " << W[0] << " " << W[1] <<  " " << W[2] << std::endl; 
    	
    	return W;
}

