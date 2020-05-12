# Astrophysical Neutrino Decay

| arXiv | Zenodo DOI |
|:-----:|:---:|
|[![arXiv](https://img.shields.io/badge/arXiv-2005.xxxxx-orange.svg)](https://arXiv.org/abs/2005.xxxxx)|[![DOI](https://zenodo.org/badge/DOI/xx.xxxx/zenodo.xxxxxxx.svg)](https://doi.org/xx.xxxx/zenodo.xxxxxxx)|

This code allows for the calculation of neutrino flavor changing probabilities in the presence of invisible or visible neutrino decay over cosmological distances. Every figure in the paper can be made via the functions in [src/Figures.cpp](src/Figures.cpp). To run one of them uncomment the relevant line in [src/main.cpp](src/main.cpp), compile, and execute `./main`.

## Usage
The main object governing the decay parameters is the struct `decay_params` found in [include/Parameters.h](include/Parameters.h) as well as `dp_bm` which is defined with the benchmark parameters. The invisible decay probability is the sum of the SM and depletion components, so to calculate the invisible decay probability composed of the SM and the depletion contributions from &nu;<sub>&mu;</sub> to &nu;<sub>e</sub> where the relevant decay parameters including the final energy are all defined in the variable `dp` of type `decay_params`, `PSM(m, e, dp) + Pdep(m, e, dp)` and for the visible case one simply addes in the regeneration term, `PSM(m, e, dp) + Pdep(m, e, dp) + Preg(m, e, dp)`. Note that the regeneration term is computationally rather slow.

There are several IceCube related functions. `Rtc(dp, vis)` calculates the track to cascade ratio for visible decay if `vis` is true and for invisible decay otherwise. The spectral indices at Earth over 100 TeV &lt; E<sub>f</sub> &lt; 1 PeV can be calculated with `IC_gamma(dp, vis, gamma_track, gamma_cascade)` which calculates the spectral index of both tracks and cascades simultaneously. Note that changing E<sub>f</sub> in `dp` does nothing here. The &chi;<sup>2</sup> test statistic can be calculated with `Chisq(dp, vis, m1min)` where the boolean `m1min` determines if the function should minimize over m<sub>1</sub> with a prior from cosmology or not. In either case, the m<sub>1</sub> prior is calculated relative to the smallest allowed value by oscillations.

## Units
- **Energy**: GeV
- **Mass**: eV
- **Distance**: km

Oscillation data is set to the best fit parameters from `nu-fit v4.1` in the normal mass ordering without Super-K data.

## Dependencies
This depends on gsl which can be found in libgsl-dev on ubuntu. This code also uses openmp for the longer calculations

## Bugs and Features
If any bugs are identified or any features suggested, please use the tools (issues, pull requests, etc.) on github.

## Reference
If you use this code please reference **[arXiv:2005.xxxxx](https://arxiv.org/abs/2005.xxxxx)** and **[doi:xx.xxxx/zenodo.xxxxxxx](https://doi.org/xx.xxxx/zenodo.xxxxxxx)**.

