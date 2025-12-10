#include <vector>
#include "lib/TransformValues.h"
#include <cmath>

std::vector<double> FindR(std::vector<double> W_0, 
                          std::vector<double> W_1, 
                          std::vector<double> W_2) {

    std::vector<double> U_0(3, 0.0), U_1(3, 0.0), U_2(3, 0.0);
    ConvertWtoU(W_0, U_0);
    ConvertWtoU(W_1, U_1);
    ConvertWtoU(W_2, U_2);

    std::vector<double> r(3, 0.0);
    for (int i = 0; i < 3; i++) {
        if (U_2[i] - U_1[i] <= 1e-12) 
            r[i] = (U_1[i] - U_0[i] <= 1e-12) ?  1 : 1e10;
        else
            r[i] = std::min(1e5, (U_1[i] - U_0[i])/(U_2[i] - U_1[i]));
    }
    

    return r;
}

double CHARM(double r) {
    
    /*if (r < 0.0) return 0.0;
    else if (r > 1e4) return 3.0;
    else return std::min(3.0, r*(3 * r + 1) / std::pow(r + 1 , 2));
    
    return std::max(0.0, std::min(1.0, r));*/

    //return (r > 1e4) ? 2 : std::max(0.0 , std::max(std::min(2*r, 1.0), std::min(r, 2.0)));
    return (r > 1e4) ? 0 : 2*r/(std::pow(r, 2.0) + 1);
    
}