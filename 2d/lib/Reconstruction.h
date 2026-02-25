#ifndef _RECONSTRUCTION_H_
#define _RECONSTRUCTION_H_

#include "Types.h"

void ReconstructGodunov(const Field& W,
                        Field& W_L,
                        Field& W_R,
                        int dir);

void ReconstructKolgan(const Field& W,
                       Field& W_L,
                       Field& W_R,
                       int dir);

void Reconstruct(const Field& W,
                 Field& W_L,
                 Field& W_R,
                 int dir);


#endif