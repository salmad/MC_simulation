#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

# include "../molecule.h"
# include "../global.h"
# include "polymer.h"
using namespace std;


// interaction potentials
inline double Distance(molecule m1, molecule m2);
inline double electrostatic_energy(molecule m1, molecule m2);
inline double lj_energy(molecule m1, molecule m2);


//  calculation of energy
double constraint_energy(molecule m1);
double calc_energy(molecule M[]);
double recalc_energy(int i, molecule M[]);

// different mc steps
int mc_step(molecule M[], int mol_ind);
double mc_steps(molecule M[],  int N);

molecule delete_molecule(molecule M[],int type);
molecule delete_molecule(molecule M[],int type,int id);

void add_molecule(molecule M[], molecule molecule_to_add);
molecule add_molecule(molecule M[], int type);

int mc_delete_ions_step(molecule M[], double mu=0.0,double N0=N_molecules/2);
double mc_delete_ions_steps(molecule M[],double mu=0.0 , double N0=N_molecules/2,  int N=1);
int mc_add_ions_step(molecule M[], double mu=0.0,double N0=N_molecules/2);
double mc_add_ions_steps(molecule M[],double mu=0.0 , double N0=N_molecules/2,  int N=1);


// analysis and plotting
int gnuplot(molecule M[],int j);
int gnuplot(polymer poly[],int j);
int radial_distribution(molecule M[], int counts[], int type);
int axis_distribution(molecule M[], int counts[], int type, int axis);
void plot_radius(int  counts[], int j);
double plot_radius(int  anion[], int cation[], int j);
void plot_axis(int  anion[], int cation[], int j);
void plot_axis(int  counts[], int j);
void return_molecules(molecule M[]);

//thermodynamic properties
double bulk_conc(double mu);
double Gamma(double mu,double N_av=(double)N_molecules);
double mc_any_steps(molecule M[],double mu=0.0,double N0=N_molecules*1.0/2.0,int N=1);
//double mc_any_steps(molecule M[],double mu=0.0,double N0=N_molecules*1.0/2.0,int N=1);
double bulk_conc2(double mu);
double Gamma2(double mu,double N_av=(double)N_molecules);




#endif // INTEGRATE_H_INCLUDED
