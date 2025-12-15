#include <vector>
#include "lib/TransformValues.h"
#include <cmath>

double beta = 1.5;

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

    if (r < 0.0) return 0.0;
    else if (r > 1e4) return 3.0;
    else return std::min(3.0, r*(3 * r + 1) / std::pow(r + 1 , 2));

}

double HCUS(double r) {
    return (r > 1e4) ? 3.0 : 1.5*(r + std::abs(r))/(r + 2);
}

double HQUICK(double r) {
    return (r > 1e4) ? 4.0 : 2*(r + std::abs(r))/(r + 3);
}

double Koren(double r) {
    return (r > 1e4) ? 2.0 : std::max(0.0, std::min(2.0, std::min(2*r, (1 + 2*r)/2)));
}

double minmodLim(double r) {
    return (r > 1e4) ? 1.0 : std::max(0.0, std::min(1.0, r));
}

double MC(double r) {
    return (r > 1e4) ? 2.0 : std::max(0.0, std::min(2.0, std::min(2*r, 0.5*(1 + r))));
}

double OsherLim(double r) {
    return (r > 1e4) ? beta : std::max(0.0, std::min(r, beta));
}

double ospre(double r) {
    return (r > 1e4) ? 1.5 : 1.5*(std::pow(r, 2) + r/(std::pow(r, 2) + r + 1));
}

double smart(double r) {
    return (r > 1e4) ? 4.0 : std::max(0.0, std::min(2*r, std::min(0.25 + 0.75*r, 4.0)));
}

double superbeeLim(double r) {
    return (r > 1e4) ? 2.0 : std::max(0.0, std::max(std::min(2*r, 1.0), std::min(r, 2.0)));
}

double Sweby(double r) {
    return (r > 1e4) ? 2.0 : std::max(0.0, std::max(std::min(beta*r, 1.0), std::min(r, beta)));
}

double UMIST(double r) {
    return (r > 1e4) ? 2.0 : std::max(0.0, std::min(std::min(2*r, 0.25 + 0.75*r), std::min(0.75 + 0.25*r, 2.0)));
}

double vanAlbada1(double r) {
    return (r > 1e4) ? 1.0 : (std::pow(r, 2) + r)/(std::pow(r, 2) + 1.0);
}

double vanAlbada2(double r) {
    return (r > 1e4) ? 0.0 : 2*r/(std::pow(r, 2) + 1.0);
}

double vanLeerLim(double r) {
    return (r > 1e4) ? 2.0 : (r + std::abs(r))/(1 + std::abs(r));
}
