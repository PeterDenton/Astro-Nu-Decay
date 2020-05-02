#ifndef Parameters_H
#define Parameters_H

#include <complex>

enum flavor {e = 0, m = 1, t = 2};
enum nunubar {nu2nu = 0, nu2nubar = 1, nu2nu_and_nu2nubar = 2};

// oscillation parameters
extern double t12, t13, t23, Dmsq21, Dmsq31, delta;
extern double Dmsq32;
extern double c12, s12, c12sq, s12sq, s212, c212;
extern double c13, s13, c13sq, s13sq, s213, c213;
extern double c23, s23, c23sq, s23sq, s223, c223;
extern std::complex<double> eid;
extern std::complex<double> U[3][3];
extern double Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;

void Recalc_Parameters();

// unit conversions
extern const double eVsqkm_to_GeV, kmGeV;

// 1502.01589
extern const double OmegaM;
extern const double OmegaLambda;
extern const double H0; // km/s/Mpc

extern const double c; // m/s
extern const double Mpc2m; // 1 pc is 3e16 m
extern const double LH; // km

// Main struct for the decay parameters
// m1: mass of lightest in eV
// gij: the off diagonal couplings, all = 0 => SM
// Ef: initial neutrino energy in GeV
// gamma: spectral index of source spectrum
// x: distance, either redshift or km as determined by x_is_redshift
// nunubar: if we are calculating the sum nu->nu + nu->nubar or if we are just calculating nu->nu
struct decay_params {double m1, g21, g31, g32, gp21, gp31, gp32, Ef, gamma, x; bool x_is_redshift; nunubar nnb;};

// benchmark
extern decay_params dp_bm;

#endif
