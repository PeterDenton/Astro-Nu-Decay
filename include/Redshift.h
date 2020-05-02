#ifndef Redshift_H
#define Redshift_H

#include "Parameters.h"

double RSE(double z);
extern double RSE_integral_ykbh, RSE_integral_hermes;
double RSE_integral(flavor a, flavor b, decay_params dp, bool ykbh);
extern const double zmax;
double h(double z);
double I2(double z);
void I2_precalc();

// same as I2 but with unit km
double Lastro(double zi, double zf);
double Lastro(double z);

#endif
