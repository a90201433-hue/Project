#ifndef _TRANSFORM_VALUES_H_
#define _TRANSFORM_VALUES_H_

void ConvertWtoU(std::vector<std::vector<double>> W,
		 std::vector<std::vector<double>>& U);

void ConvertUtoW(std::vector<std::vector<double>>& W,
		 std::vector<std::vector<double>> U);

#endif
