#include <vector>
#include <cmath>
#include <iostream>

extern int Nx, fict;


double L1(std::vector<double> U) {
    double Norm = 0;
    for(int i = fict; i < Nx + fict; i++)
    {
        Norm += std::abs(U[i]);
    }

    return Norm;
}

double LInf(std::vector<double> U) {
    double Norm = 0;
    for(int i = fict; i < Nx + fict; i++)
    {
        if (Norm < std::abs(U[i]))
        {
            Norm = std::abs(U[i]);
        }
    }
    return Norm;
}

double RatioL(std::vector<double> U) {
    double c;
    c = L1(U)/LInf(U);
    return c;
}