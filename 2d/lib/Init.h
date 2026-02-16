#ifndef _INIT_H_
#define _INIT_H_

#include "Types.h"

void readConfig(const std::string& config_path);

void Grid(std::vector<double>& x, std::vector<double>& xc,
		  std::vector<double>& y, std::vector<double>& yc);

void InitialZeros(std::vector<std::vector<double>> &W_zeros, int intern_size);

void InitValues(Field& W, 
				const std::vector<double>& x, 
				const std::vector<double>& y,
				const std::string& config_path);

#endif
