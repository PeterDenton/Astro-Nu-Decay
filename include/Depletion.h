#ifndef Depletion_H
#define Depletion_H

#include "Parameters.h"

double Gij(int i, int j, double Ei, decay_params dp);

double Pdep(flavor a, flavor b, decay_params dp);

#endif
