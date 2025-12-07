#ifndef _TRANSFORM_VALUES_H_
#define _TRANSFORM_VALUES_H_

void ConvertWtoU(std::vector<std::vector<double>> W,
		 std::vector<std::vector<double>>& U);

void ConvertWtoU(std::vector<double> W,
		 std::vector<double>& U);

void ConvertUtoW(std::vector<std::vector<double>>& W,
		 std::vector<std::vector<double>> U);

void ConvertUtoW(std::vector<double>& W,
		 std::vector<double> U);

void FindFfromW(std::vector<double>& F, 
		std::vector<double> W);

#endif
