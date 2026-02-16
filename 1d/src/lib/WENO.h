#ifndef _WENO_H
#define _WENO_H

using DataArray = std::vector<std::vector<double>>;

std::vector<double> WENOWeight(double f1, double f2, double f3, double f4, double f5);

void WENO(
	DataArray U,
	DataArray& U_L,
	DataArray& U_R);

void WENOStreams(
	DataArray W,
	DataArray& F);

#endif
