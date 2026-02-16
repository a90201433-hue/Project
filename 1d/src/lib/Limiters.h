#ifndef _LIMITERS_H_
#define _LIMITERS_H_

#include <functional>

using LimiterFunction = std::function<double(double)>;

std::vector<double> FindR(std::vector<double> W_0, 
                          std::vector<double> W_1, 
                          std::vector<double> W_2);

double CHARM(double r);

double HCUS(double r);

double HQUICK(double r);

double Koren(double r);

double minmodLim(double r);

double MC(double r);

double OsherLim(double r);

double ospre(double r);

double smart(double r);

double superbeeLim(double r);

double Sweby(double r);

double UMIST(double r);

double vanAlbada1(double r);

double vanAlbada2(double r);

double vanLeerLim(double r);

#endif
