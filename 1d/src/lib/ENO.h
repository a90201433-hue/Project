#ifndef _ENO_H_
#define _ENO_H

using DataArray = std::vector<std::vector<double>>;

double ENOPolinom(double f1, double f2, double f3, double x);

void ENO(
	DataArray U,
	DataArray Slope,
	DataArray& U_L,
	DataArray& U_R);

#endif
