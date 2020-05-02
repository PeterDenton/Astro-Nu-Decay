#include "IceCube.h"
#include "Parameters.h"
#include "SM.h"
#include "Depletion.h"
#include "Regeneration.h"
#include "Power_Law_Fit.h"

// Track to cascade ratio
double Rtc(decay_params dp, bool visible)
{
	// Calculate the regeneration part if necessary
	double P21, P31, P32;
	if (visible)
	{
		P21 = Pregij(2, 1, dp);
		P31 = Pregij(3, 1, dp);
		P32 = Pregij(3, 2, dp);
	}
	else
	{
		P21 = 0;
		P31 = 0;
		P32 = 0;
	}

	// Calculate the six relevant probabilities
	double Pee, Pem, Pet, Pme, Pmm, Pmt;
	Pee = PSM(e, e, dp) + Pdep(e, e, dp) +
			+ std::norm(U[e][1]) * std::norm(U[e][0]) * P21
			+ std::norm(U[e][2]) * std::norm(U[e][0]) * P31
			+ std::norm(U[e][2]) * std::norm(U[e][1]) * P32;
	Pem = PSM(e, m, dp) + Pdep(e, m, dp) +
			+ std::norm(U[e][1]) * std::norm(U[m][0]) * P21
			+ std::norm(U[e][2]) * std::norm(U[m][0]) * P31
			+ std::norm(U[e][2]) * std::norm(U[m][1]) * P32;
	Pet = PSM(e, t, dp) + Pdep(e, t, dp) +
			+ std::norm(U[e][1]) * std::norm(U[t][0]) * P21
			+ std::norm(U[e][2]) * std::norm(U[t][0]) * P31
			+ std::norm(U[e][2]) * std::norm(U[t][1]) * P32;
	Pme = PSM(m, e, dp) + Pdep(m, e, dp) +
			+ std::norm(U[m][1]) * std::norm(U[e][0]) * P21
			+ std::norm(U[m][2]) * std::norm(U[e][0]) * P31
			+ std::norm(U[m][2]) * std::norm(U[e][1]) * P32;
	Pmm = PSM(m, m, dp) + Pdep(m, m, dp) +
			+ std::norm(U[m][1]) * std::norm(U[m][0]) * P21
			+ std::norm(U[m][2]) * std::norm(U[m][0]) * P31
			+ std::norm(U[m][2]) * std::norm(U[m][1]) * P32;
	Pmt = PSM(m, t, dp) + Pdep(m, t, dp) +
			+ std::norm(U[m][1]) * std::norm(U[t][0]) * P21
			+ std::norm(U[m][2]) * std::norm(U[t][0]) * P31
			+ std::norm(U[m][2]) * std::norm(U[t][1]) * P32;

    double Nt, Nc;
	// Pion decay
	Nt = Pem + 2 * Pmm;
	Nc = Pee + 2 * Pme + Pet + 2 * Pmt;

    return Nt / Nc;
}

// Calculate the spectral index at the Earth over the range 100 TeV to 1 PeV assuming pion decay sources
void IC_gamma(decay_params dp, bool visible, double &gamma_track, double &gamma_cascade)
{
	double Efmin, Efmax, Efscale;
	const int n = 20; // n=10 => <0.5% error, n=20 => <0.3% error
	double Efs[n + 1], Phi_tracks[n + 1], Phi_cascades[n + 1];

	// The range that IceCube is sensitive to
	Efmin = 1e5;
	Efmax = 1e6;
	Efscale = pow(Efmax / Efmin, 1. / n);

	double p_track, p_cascade, p21, p31, p32, pee, pem, pet, pme, pmm, pmt;
	for (int i = 0; i <= n; i++)
	{
		dp.Ef = Efmin * pow(Efscale, i);
		Efs[i] = dp.Ef;

		// Track probability
		p_track = 0;
		p_track += PSM(e, m, dp) + Pdep(e, m, dp); // e->m
		p_track += 2 * (PSM(m, m, dp) + Pdep(m, m, dp)); // m->m

		// Cascade probability
		p_cascade = 0;
		p_cascade += PSM(e, e, dp) + Pdep(e, e, dp); // e->e
		p_cascade += PSM(e, t, dp) + Pdep(e, t, dp); // e->t
		p_cascade += 2 * (PSM(m, e, dp) + Pdep(m, e, dp)); // m->e
		p_cascade += 2 * (PSM(m, t, dp) + Pdep(m, t, dp)); // m->t

		// Include the regeneration terms if necessary
		if (visible)
		{
			// Calculate the three mass basis probabilities
			p21 = Pregij(2, 1, dp);
			p31 = Pregij(3, 1, dp);
			p32 = Pregij(3, 2, dp);

			// Calculate the six relevant flavor basis probabilities
			pee = Ue2sq * Ue1sq * p21 + Ue3sq * Ue1sq * p31 + Ue3sq * Ue2sq * p32;
			pem = Ue2sq * Um1sq * p21 + Ue3sq * Um1sq * p31 + Ue3sq * Um2sq * p32;
			pet = Ue2sq * Ut1sq * p21 + Ue3sq * Ut1sq * p31 + Ue3sq * Ut2sq * p32;
			pme = Um2sq * Ue1sq * p21 + Um3sq * Ue1sq * p31 + Um3sq * Ue2sq * p32;
			pmm = Um2sq * Um1sq * p21 + Um3sq * Um1sq * p31 + Um3sq * Um2sq * p32;
			pmt = Um2sq * Ut1sq * p21 + Um3sq * Ut1sq * p31 + Um3sq * Ut2sq * p32;

			p_track += pem + 2 * pmm; // track: e->m and m->m
			p_cascade += pee + pet + 2 * (pme + pmt); // cascade: e->e, e->t, m->e, m->t
		}

		// Multiply back by the initial spectra to get the observed spectrum
		Phi_tracks[i] = p_track * Spectrum(dp.Ef, dp.gamma, 0);
		Phi_cascades[i] = p_cascade * Spectrum(dp.Ef, dp.gamma, 0);
	} // i, n, Ef

	// Call the power law fitting function to get the exponent, discard the normalization (tmp)
	double tmp;
	Power_Law_Fit(Efs, Phi_tracks, n, 1e5, tmp, gamma_track);
	Power_Law_Fit(Efs, Phi_cascades, n, 1e5, tmp, gamma_cascade);
}

