#ifndef _MADER_SOLVER_H_
#define _MADER_SOLVER_H_

#include "Types.h"

extern Field omega;

void InitializeMaderSolver(int Nx,int Ny);

void MaderTimeStep(
Field& W,
Field& W_new,
double dt,
double dx,
double dy
);

#endif