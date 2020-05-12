#ifndef IceCube_H
#define IceCube_H

#include "Parameters.h"

double Rtc(decay_params dp, bool visible = true);
double IC_gamma(flavor a, flavor b, decay_params dp, bool visible = true);
void IC_gamma(decay_params dp, bool visible, double &gamma_track, double &gamma_cascade);
double Chisq(decay_params dp, bool visible, bool m1min);

#endif
