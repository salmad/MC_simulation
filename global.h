#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED
# include <string>
#include <stdio.h>
#include <stdlib.h>
# include <cmath>

using namespace std;

extern string foldername;
//time steps and number of bins
extern const double time_step;
extern const int nbins;
//Droplet size
extern const int geometry;  // 0 if cubic , 1 if spherical
extern const double L; // Nanometers
extern const double V;
extern const double R;


//Electrostatic parameter
extern const double bjerrum;

//Lennard-jones parameters
extern const double lj_eps;
extern const double lj_sigma;
extern const double lj_shift;
extern const double lj_cut;
extern const double r_off_anions;

extern const int zero_adsorbption;
extern const double lj_eps_constraint_anions;
extern const double r_off_anions;
//const double lj_shift_anions=0.25;
//const double lj_cut_constraint_anions=1.12246204831*lj_sigma;
extern const double lj_shift_anions;
extern const double lj_cut_constraint_anions;
extern const double lj_shift_cations;

extern int N_molecules;
extern int N_polymers;
extern int mol_id;
extern int pol_id;
extern const double kT;
extern const double lambda;
extern double k;



#endif // GLOBAL_H_INCLUDED
