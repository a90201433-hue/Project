#ifndef _ACOUSTIC_H_
#define _ACOUSTIC_H_

void AcoustPStarAndUStar(std::vector<double> W_L, 
			 std::vector<double> W_R, 
			 std::vector<double>& W_ac_star);

std::vector<double> ConfigurationAcousticWaves(
		std::vector<double> W_L, 
		std::vector<double> W_R, 
		std::vector<double> W_ac_star, 
		
		double x, double t);

#endif
