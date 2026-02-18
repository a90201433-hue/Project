#ifndef _FLIC_H_
#define _FLIC_H

#include <functional>
#include "Types.h"

void FLIC_L(Field W,
            Field& W_tilde,
            std::vector<double> x,
            std::vector<double> y,
            double dt,
            int dir);

void FLIC_E(Field W_tilde,
            Field& W_new,
            std::vector<double> x,
            std::vector<double> y,
            double dt,
            int dir);

void FLIC(Field& W_new, 
		  const Field& W,
          const std::vector<double>& x, 
		  const std::vector<double>& y, 
		  double dt);

#endif