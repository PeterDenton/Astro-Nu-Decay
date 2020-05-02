#ifndef Regeneration_H
#define Regeneration_H

#include "Parameters.h"

double Spectrum(double Ei, double gamma, double z);
double Pregij(int i, int j, decay_params dp);
double Preg(flavor a, flavor b, decay_params dp);

#endif
